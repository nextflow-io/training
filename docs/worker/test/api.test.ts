import { SignJWT, exportJWK, generateKeyPair } from 'jose';
import { beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import worker, { type Env } from '../src/index';
import { returnPath } from '../src/api';

const ISSUER = 'https://login.example.com/';
const rows = new Map<string, string>();
const db = {
  prepare: (sql: string) => ({
    bind: (...args: string[]) => ({
      first: async () => (sql.startsWith('SELECT') ? (rows.get(args[0]) ?? null) : null),
      run: async () => void rows.set(args[0], args[1]),
    }),
  }),
};
const prod = {
  BUCKET: {},
  DB: db,
  OIDC_ISSUER: ISSUER,
  OIDC_CLIENT_ID: 'client',
  OIDC_CLIENT_SECRET: 'shh',
  AUTH_SECRET: 'x'.repeat(32),
} as unknown as Env;
const ctx = { waitUntil() {} } as unknown as ExecutionContext;

let signer: CryptoKey;
let idTokenClaims: Record<string, unknown>;
let cacheUsed: boolean;

beforeAll(async () => {
  const { privateKey, publicKey } = await generateKeyPair('RS256');
  signer = privateKey;
  const jwk = { ...(await exportJWK(publicKey)), kid: 'k1', alg: 'RS256' };
  vi.stubGlobal('fetch', async (input: RequestInfo | URL, init?: RequestInit) => {
    const url = String(input instanceof Request ? input.url : input);
    if (url === `${ISSUER}.well-known/openid-configuration`) {
      return Response.json({
        issuer: ISSUER,
        authorization_endpoint: `${ISSUER}authorize`,
        token_endpoint: `${ISSUER}oauth/token`,
        jwks_uri: `${ISSUER}.well-known/jwks.json`,
      });
    }
    if (url === `${ISSUER}.well-known/jwks.json`) return Response.json({ keys: [jwk] });
    if (url === `${ISSUER}oauth/token`) {
      const body = new URLSearchParams(String(init?.body));
      if (body.get('code') !== 'good' || body.get('code_verifier') !== 'v') return new Response('bad', { status: 400 });
      const { aud = 'client', ...claims } = idTokenClaims;
      const idToken = await new SignJWT(claims)
        .setProtectedHeader({ alg: 'RS256', kid: 'k1' })
        .setIssuer(ISSUER)
        .setAudience(aud as string)
        .setSubject('seqera|123')
        .setIssuedAt()
        .setExpirationTime('5m')
        .sign(signer);
      return Response.json({ id_token: idToken });
    }
    return new Response('unexpected fetch', { status: 500 });
  });
});

beforeEach(() => {
  rows.clear();
  idTokenClaims = { name: 'Ada Lovelace', email: 'ada@example.com' };
  cacheUsed = false;
  vi.stubGlobal('caches', {
    default: {
      match: async () => {
        cacheUsed = true;
      },
      put: async () => {
        cacheUsed = true;
      },
    },
  });
});

const call = (
  path: string,
  init: { method?: string; body?: string; headers?: Record<string, string> } = {},
  env = prod,
  origin = 'https://training.nextflow.io',
) => worker.fetch(new Request(`${origin}${path}`, init), env, ctx);

const cookies = (res: Response) =>
  Object.fromEntries(res.headers.getSetCookie().map((c) => c.split(';')[0].split(/=(.*)/s).slice(0, 2)));

async function logIn(): Promise<string> {
  const start = await call('/api/auth/login?return=/latest/hello_nextflow/');
  const oauth = cookies(start).nf_oauth;
  const state = new URL(start.headers.get('Location')!).searchParams.get('state');
  const saved = JSON.parse(decodeURIComponent(oauth));
  const cookie = `nf_oauth=${encodeURIComponent(JSON.stringify({ ...saved, verifier: 'v' }))}`;
  const res = await call(`/api/auth/callback?code=good&state=${state}`, { headers: { Cookie: cookie } });
  expect(res.headers.get('Location')).toBe('https://training.nextflow.io/latest/hello_nextflow/');
  return `nf_session=${cookies(res).nf_session}`;
}

describe('auth', () => {
  it('starts the code flow with PKCE', async () => {
    const res = await call('/api/auth/login?return=/latest/');
    expect(res.status).toBe(302);
    const to = new URL(res.headers.get('Location')!);
    expect(to.origin + to.pathname).toBe(`${ISSUER}authorize`);
    expect(to.searchParams.get('code_challenge_method')).toBe('S256');
    expect(to.searchParams.get('redirect_uri')).toBe('https://training.nextflow.io/api/auth/callback');
    expect(to.searchParams.get('scope')).toBe('openid profile email');
    expect(res.headers.getSetCookie()[0]).toMatch(
      /^nf_oauth=.*Path=\/api\/auth; Max-Age=600; SameSite=Lax; HttpOnly; Secure$/,
    );
  });

  it('logs in, setting the session and the page marker', async () => {
    const start = await call('/api/auth/login');
    const state = new URL(start.headers.get('Location')!).searchParams.get('state');
    const saved = JSON.parse(decodeURIComponent(cookies(start).nf_oauth));
    const cookie = `nf_oauth=${encodeURIComponent(JSON.stringify({ ...saved, verifier: 'v' }))}`;
    const res = await call(`/api/auth/callback?code=good&state=${state}`, { headers: { Cookie: cookie } });
    expect(res.status).toBe(302);
    expect(res.headers.get('Location')).toBe('https://training.nextflow.io/latest/');
    const set = res.headers.getSetCookie();
    expect(set.find((c) => c.startsWith('nf_session='))).toMatch(
      /Path=\/api; Max-Age=2592000; SameSite=Lax; HttpOnly; Secure/,
    );
    expect(set.find((c) => c.startsWith('nf_login=1'))).not.toMatch(/HttpOnly/);
    expect(set.find((c) => c.startsWith('nf_oauth='))).toMatch(/Max-Age=0/);
  });

  it('rejects a callback whose state does not match', async () => {
    const cookie = `nf_oauth=${encodeURIComponent(JSON.stringify({ state: 'a', verifier: 'v', back: '/' }))}`;
    expect((await call('/api/auth/callback?code=good&state=b', { headers: { Cookie: cookie } })).status).toBe(400);
    expect((await call('/api/auth/callback?code=good&state=a')).status).toBe(400);
  });

  it('falls back to preferred_username, then email, then sub for the display name', async () => {
    const name = async (claims: Record<string, unknown>) => {
      idTokenClaims = claims;
      const cookie = await logIn();
      const res = await call('/api/progress', { headers: { Cookie: cookie } });
      return ((await res.json()) as { user: { name: string } }).user.name;
    };
    expect(await name({ preferred_username: 'ada', email: 'ada@example.com' })).toBe('ada');
    expect(await name({ email: 'ada@example.com' })).toBe('ada@example.com');
    expect(await name({})).toBe('seqera|123');
  });

  it('rejects an ID token for another client', async () => {
    idTokenClaims = { ...idTokenClaims, aud: 'someone-else' };
    const cookie = `nf_oauth=${encodeURIComponent(JSON.stringify({ state: 'a', verifier: 'v', back: '/' }))}`;
    expect((await call('/api/auth/callback?code=good&state=a', { headers: { Cookie: cookie } })).status).toBe(401);
  });

  it('only returns to paths on this site', () => {
    expect(returnPath('/latest/x/?a=1')).toBe('/latest/x/?a=1');
    expect(returnPath('//evil.com/')).toBe('/latest/');
    expect(returnPath('/\\evil.com/')).toBe('/latest/');
    expect(returnPath('https://evil.com/')).toBe('/latest/');
    expect(returnPath(null)).toBe('/latest/');
  });

  it('logs out', async () => {
    const res = await call('/api/auth/logout?return=//evil.com/');
    expect(res.headers.get('Location')).toBe('https://training.nextflow.io/latest/');
    expect(res.headers.getSetCookie().every((c) => c.includes('Max-Age=0'))).toBe(true);
  });

  it('signs in a dev user on localhost only, when no issuer is set', async () => {
    const dev = { BUCKET: {}, DB: db } as unknown as Env;
    const local = await call('/api/auth/login?return=/x/', {}, dev, 'http://localhost:4321');
    expect(local.headers.get('Location')).toBe('http://localhost:4321/x/');
    expect(local.headers.getSetCookie().join()).not.toMatch(/Secure/);
    expect((await call('/api/auth/login', {}, dev)).status).toBe(503);
    expect((await call('/api/auth/login', {}, { ...dev, AUTH_SECRET: 'x'.repeat(32) } as Env)).status).toBe(503);
    expect((await call('/api/auth/login', {}, prod, 'http://localhost:4321')).headers.get('Location')).toMatch(
      /^https:\/\/login\.example\.com\//,
    );
  });
});

describe('progress', () => {
  const put = (body: string, cookie: string, headers: Record<string, string> = {}) =>
    call('/api/progress', {
      method: 'PUT',
      body,
      headers: {
        Cookie: cookie,
        Origin: 'https://training.nextflow.io',
        'Content-Type': 'application/json',
        ...headers,
      },
    });

  it('needs a session', async () => {
    const res = await call('/api/progress');
    expect(res.status).toBe(401);
    expect(res.headers.getSetCookie()[0]).toMatch(/^nf_login=; .*Max-Age=0/);
    expect((await call('/api/progress', { headers: { Cookie: 'nf_session=forged' } })).status).toBe(401);
    expect((await put('{}', '')).status).toBe(401);
  });

  it('saves and loads the state', async () => {
    const cookie = await logIn();
    const empty = await call('/api/progress', { headers: { Cookie: cookie } });
    expect(await empty.json()).toEqual({
      user: { sub: 'seqera|123', name: 'Ada Lovelace', email: 'ada@example.com' },
      state: null,
    });
    expect(empty.headers.get('Cache-Control')).toBe('no-store');
    expect((await put('{"basket":["hello_nextflow"]}', cookie)).status).toBe(204);
    const res = await call('/api/progress', { headers: { Cookie: cookie } });
    expect(((await res.json()) as { state: unknown }).state).toEqual({ basket: ['hello_nextflow'] });
    expect([...rows.keys()]).toEqual(['seqera|123']);
  });

  it('validates saves', async () => {
    const cookie = await logIn();
    expect((await put('{}', cookie, { Origin: 'https://evil.com' })).status).toBe(403);
    expect((await put('{}', cookie, { Origin: '' })).status).toBe(403);
    expect((await put('{}', cookie, { 'Content-Type': 'text/plain' })).status).toBe(415);
    expect((await put(`{"x":"${'a'.repeat(40000)}"}`, cookie)).status).toBe(413);
    expect((await put('[1]', cookie)).status).toBe(400);
    expect((await put('null', cookie)).status).toBe(400);
    expect((await put('{', cookie)).status).toBe(400);
    expect(rows.size).toBe(0);
  });

  it('bypasses the cache', async () => {
    const cookie = await logIn();
    await call('/api/progress', { headers: { Cookie: cookie } });
    expect(cacheUsed).toBe(false);
  });

  it('is unavailable without a session secret', async () => {
    expect((await call('/api/progress', {}, { ...prod, AUTH_SECRET: undefined } as Env)).status).toBe(503);
  });
});
