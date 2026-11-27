// Optional Seqera login (OIDC code flow with PKCE) and synced My courses progress.
import { SignJWT, createRemoteJWKSet, jwtVerify } from 'jose';

export interface ApiEnv {
  DB: D1Database;
  OIDC_ISSUER?: string;
  OIDC_CLIENT_ID?: string;
  OIDC_CLIENT_SECRET?: string;
  AUTH_SECRET?: string;
}

type User = { sub: string; name: string; email: string };

const SESSION = 'nf_session';
/** Readable by page scripts: tells them to sync, so logged-out visitors make no API calls. */
const MARKER = 'nf_login';
const OAUTH = 'nf_oauth';
const MAX_AGE = 60 * 60 * 24 * 30;
const MAX_BODY = 32 * 1024;
const DEV_USER: User = { sub: 'dev|learner', name: 'Dev Learner', email: 'dev@localhost' };

export async function api(request: Request, env: ApiEnv): Promise<Response> {
  const url = new URL(request.url);
  const route = `${request.method} ${url.pathname.replace(/\/$/, '')}`;
  const dev = !env.OIDC_ISSUER && ['localhost', '127.0.0.1'].includes(url.hostname);
  const secret = env.AUTH_SECRET ?? (dev ? 'local-development-only-secret' : undefined);
  if (!secret) return text('Login is not configured', 503);
  const key = new TextEncoder().encode(secret);

  switch (route) {
    case 'GET /api/auth/login': {
      const back = returnPath(url.searchParams.get('return'));
      if (dev) return signIn(url, DEV_USER, key, back);
      if (!env.OIDC_ISSUER || !env.OIDC_CLIENT_ID) return text('Login is not configured', 503);
      const oidc = await discover(env.OIDC_ISSUER);
      const state = random();
      const verifier = random();
      const target = new URL(oidc.authorization_endpoint);
      target.search = new URLSearchParams({
        response_type: 'code',
        client_id: env.OIDC_CLIENT_ID,
        redirect_uri: `${url.origin}/api/auth/callback`,
        scope: 'openid profile email',
        state,
        code_challenge: base64url(await crypto.subtle.digest('SHA-256', new TextEncoder().encode(verifier))),
        code_challenge_method: 'S256',
      }).toString();
      const res = redirect(target.href);
      res.headers.append(
        'Set-Cookie',
        cookie(url, OAUTH, encodeURIComponent(JSON.stringify({ state, verifier, back })), 600, '/api/auth'),
      );
      return res;
    }

    case 'GET /api/auth/callback': {
      let saved: { state?: string; verifier?: string; back?: string } = {};
      try {
        saved = JSON.parse(decodeURIComponent(readCookie(request, OAUTH) ?? '{}'));
      } catch {}
      const code = url.searchParams.get('code');
      if (!code || !saved.state || !saved.verifier || url.searchParams.get('state') !== saved.state) {
        return text('Login failed: the request has expired or did not start here. Try again.', 400);
      }
      if (!env.OIDC_ISSUER || !env.OIDC_CLIENT_ID || !env.OIDC_CLIENT_SECRET) {
        return text('Login is not configured', 503);
      }
      const oidc = await discover(env.OIDC_ISSUER);
      const tokens = await fetch(oidc.token_endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/x-www-form-urlencoded', Accept: 'application/json' },
        body: new URLSearchParams({
          grant_type: 'authorization_code',
          code,
          code_verifier: saved.verifier,
          redirect_uri: `${url.origin}/api/auth/callback`,
          client_id: env.OIDC_CLIENT_ID,
          client_secret: env.OIDC_CLIENT_SECRET,
        }),
      });
      const idToken = tokens.ok ? ((await tokens.json()) as { id_token?: string }).id_token : undefined;
      if (!idToken) return text('Login failed', 502);
      let user: User;
      try {
        const { payload } = await jwtVerify(idToken, jwks(oidc.jwks_uri), {
          issuer: oidc.issuer,
          audience: env.OIDC_CLIENT_ID,
        });
        if (!payload.sub) throw new Error('No sub');
        const str = (v: unknown) => (typeof v === 'string' ? v : '');
        const email = str(payload.email);
        const name = str(payload.name) || str(payload.preferred_username) || email || payload.sub;
        user = { sub: payload.sub, name, email };
      } catch {
        return text('Login failed', 401);
      }
      const res = await signIn(url, user, key, returnPath(saved.back ?? null));
      res.headers.append('Set-Cookie', cookie(url, OAUTH, '', 0, '/api/auth'));
      return res;
    }

    case 'GET /api/auth/logout': {
      const res = redirect(new URL(returnPath(url.searchParams.get('return')), url).href);
      res.headers.append('Set-Cookie', cookie(url, SESSION, '', 0, '/api'));
      res.headers.append('Set-Cookie', cookie(url, MARKER, '', 0, '/', false));
      return res;
    }

    case 'GET /api/progress': {
      const user = await session(request, key);
      if (!user) return loggedOut(url);
      const state = await env.DB.prepare('SELECT state FROM progress WHERE sub = ?1')
        .bind(user.sub)
        .first<string>('state');
      return json({ user, state: state ? JSON.parse(state) : null });
    }

    case 'PUT /api/progress': {
      const user = await session(request, key);
      if (!user) return loggedOut(url);
      if (request.headers.get('Origin') !== url.origin) return text('Forbidden', 403);
      if (!request.headers.get('Content-Type')?.startsWith('application/json')) return text('Expected JSON', 415);
      if (Number(request.headers.get('Content-Length') ?? 0) > MAX_BODY) return text('Too large', 413);
      const body = await request.text();
      if (body.length > MAX_BODY) return text('Too large', 413);
      let state: unknown;
      try {
        state = JSON.parse(body);
      } catch {}
      if (!state || typeof state !== 'object' || Array.isArray(state)) return text('Expected a JSON object', 400);
      await env.DB.prepare(
        'INSERT INTO progress (sub, state, updated_at) VALUES (?1, ?2, ?3) ' +
          'ON CONFLICT (sub) DO UPDATE SET state = excluded.state, updated_at = excluded.updated_at',
      )
        .bind(user.sub, body, new Date().toISOString())
        .run();
      return new Response(null, { status: 204, headers: { 'Cache-Control': 'no-store' } });
    }
  }
  return text('Not Found', 404);
}

/** Only paths on this site, so the login can't be used to send people elsewhere. */
export function returnPath(value: string | null): string {
  return value && value.startsWith('/') && !value.startsWith('//') && !value.startsWith('/\\') ? value : '/latest/';
}

async function signIn(url: URL, user: User, key: Uint8Array, back: string): Promise<Response> {
  const token = await new SignJWT({ name: user.name, email: user.email })
    .setProtectedHeader({ alg: 'HS256' })
    .setSubject(user.sub)
    .setIssuer('nextflow-training')
    .setIssuedAt()
    .setExpirationTime(`${MAX_AGE}s`)
    .sign(key);
  const res = redirect(new URL(back, url).href);
  res.headers.append('Set-Cookie', cookie(url, SESSION, token, MAX_AGE, '/api'));
  res.headers.append('Set-Cookie', cookie(url, MARKER, '1', MAX_AGE, '/', false));
  return res;
}

async function session(request: Request, key: Uint8Array): Promise<User | null> {
  const token = readCookie(request, SESSION);
  if (!token) return null;
  try {
    const { payload } = await jwtVerify(token, key, { issuer: 'nextflow-training', algorithms: ['HS256'] });
    return { sub: payload.sub!, name: String(payload.name ?? ''), email: String(payload.email ?? '') };
  } catch {
    return null;
  }
}

/** Also clears the marker, so a page stops asking once the session has expired. */
function loggedOut(url: URL): Response {
  const res = text('Not logged in', 401);
  res.headers.append('Set-Cookie', cookie(url, MARKER, '', 0, '/', false));
  return res;
}

type Discovery = { issuer: string; authorization_endpoint: string; token_endpoint: string; jwks_uri: string };
const discovered = new Map<string, Promise<Discovery>>();

function discover(issuer: string): Promise<Discovery> {
  const base = issuer.replace(/\/?$/, '/');
  let found = discovered.get(base);
  if (!found) {
    found = fetch(`${base}.well-known/openid-configuration`).then((r) => {
      if (!r.ok) throw new Error(`OIDC discovery failed: ${r.status}`);
      return r.json() as Promise<Discovery>;
    });
    found.catch(() => discovered.delete(base));
    discovered.set(base, found);
  }
  return found;
}

const keySets = new Map<string, ReturnType<typeof createRemoteJWKSet>>();
function jwks(uri: string) {
  let set = keySets.get(uri);
  if (!set) keySets.set(uri, (set = createRemoteJWKSet(new URL(uri))));
  return set;
}

function cookie(url: URL, name: string, value: string, maxAge: number, path: string, httpOnly = true): string {
  return [
    `${name}=${value}`,
    `Path=${path}`,
    `Max-Age=${maxAge}`,
    'SameSite=Lax',
    httpOnly && 'HttpOnly',
    url.protocol === 'https:' && 'Secure',
  ]
    .filter(Boolean)
    .join('; ');
}

function readCookie(request: Request, name: string): string | undefined {
  for (const part of (request.headers.get('Cookie') ?? '').split(';')) {
    const [k, ...v] = part.trim().split('=');
    if (k === name) return v.join('=');
  }
}

function random(): string {
  return base64url(crypto.getRandomValues(new Uint8Array(32)));
}

function base64url(bytes: ArrayBuffer | Uint8Array): string {
  return btoa(String.fromCharCode(...new Uint8Array(bytes)))
    .replace(/\+/g, '-')
    .replace(/\//g, '_')
    .replace(/=+$/, '');
}

function redirect(location: string): Response {
  return new Response(null, { status: 302, headers: { Location: location, 'Cache-Control': 'no-store' } });
}

function json(body: unknown): Response {
  return Response.json(body, { headers: { 'Cache-Control': 'no-store' } });
}

function text(body: string, status: number): Response {
  return new Response(body, {
    status,
    headers: { 'Content-Type': 'text/plain; charset=utf-8', 'Cache-Control': 'no-store' },
  });
}
