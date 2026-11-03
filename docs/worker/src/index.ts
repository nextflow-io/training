export interface Env {
  BUCKET: R2Bucket;
}

const PREFIX = 'training';

// Root paths from before versioning (gh-pages/_redirects).
const LEGACY = new Set([
  'advanced',
  'assets',
  'basic_training',
  'envsetup',
  'es',
  'fr',
  'hands_on',
  'hello_nextflow',
  'help',
  'nf_customize',
  'nf_develop',
  'pt',
  'troubleshoot',
]);

const TYPES: Record<string, string> = {
  html: 'text/html; charset=utf-8',
  css: 'text/css; charset=utf-8',
  js: 'text/javascript; charset=utf-8',
  mjs: 'text/javascript; charset=utf-8',
  json: 'application/json; charset=utf-8',
  xml: 'application/xml; charset=utf-8',
  txt: 'text/plain; charset=utf-8',
  md: 'text/markdown; charset=utf-8',
  svg: 'image/svg+xml',
  png: 'image/png',
  jpg: 'image/jpeg',
  jpeg: 'image/jpeg',
  gif: 'image/gif',
  webp: 'image/webp',
  avif: 'image/avif',
  ico: 'image/x-icon',
  woff: 'font/woff',
  woff2: 'font/woff2',
  ttf: 'font/ttf',
  pdf: 'application/pdf',
  wasm: 'application/wasm',
  gz: 'application/gzip',
  zip: 'application/zip',
};

export default {
  async fetch(request, env, ctx) {
    if (request.method !== 'GET' && request.method !== 'HEAD') {
      return new Response('Method Not Allowed', { status: 405, headers: { Allow: 'GET, HEAD' } });
    }
    const url = new URL(request.url);
    const key = new Request(url);
    let response = await caches.default.match(key);
    if (!response) {
      response = await route(url, env);
      if (response.status === 200) ctx.waitUntil(caches.default.put(key, response.clone()));
    }
    return request.method === 'HEAD' ? new Response(null, response) : response;
  },
} satisfies ExportedHandler<Env>;

async function route(url: URL, env: Env): Promise<Response> {
  const path = url.pathname;
  const first = path.split('/')[1];
  if (path === '/') return redirect(url, '/latest/', 302);
  if (LEGACY.has(first)) return redirect(url, `/latest${path}`, 301);

  let key: string;
  try {
    key = PREFIX + decodeURIComponent(path);
  } catch {
    return new Response('Bad Request', { status: 400 });
  }
  const object = await env.BUCKET.get(path.endsWith('/') ? `${key}index.html` : key);
  if (object) return serve(object, path, 200);
  if (!path.endsWith('/') && (await env.BUCKET.head(`${key}/index.html`))) return redirect(url, `${path}/`, 301);

  const notFound =
    (await env.BUCKET.get(`${PREFIX}/${first}/404.html`)) ?? (await env.BUCKET.get(`${PREFIX}/latest/404.html`));
  if (notFound) return serve(notFound, '/404.html', 404);
  return new Response('Not Found', { status: 404, headers: { 'X-Robots-Tag': 'noindex' } });
}

function serve(object: R2ObjectBody, path: string, status: number): Response {
  const headers = new Headers();
  const file = path.endsWith('/') ? 'index.html' : path.slice(path.lastIndexOf('/') + 1);
  const ext = file.includes('.') ? file.slice(file.lastIndexOf('.') + 1).toLowerCase() : '';
  headers.set('Content-Type', TYPES[ext] ?? object.httpMetadata?.contentType ?? 'application/octet-stream');
  headers.set('ETag', object.httpEtag);
  if (!path.startsWith('/latest/')) headers.set('X-Robots-Tag', 'noindex');
  if (status === 200) headers.set('Cache-Control', cacheControl(path));
  return new Response(object.body, { status, headers });
}

// No purge step: alias prefixes are rewritten on every deploy, so they expire quickly.
function cacheControl(path: string): string {
  if (path.includes('/_astro/')) return 'public, max-age=31536000, immutable';
  return /^\/\d/.test(path) ? 'public, max-age=86400' : 'public, max-age=300';
}

function redirect(url: URL, path: string, status: number): Response {
  return new Response(null, { status, headers: { Location: new URL(path + url.search, url).href } });
}
