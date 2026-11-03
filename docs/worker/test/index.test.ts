import { beforeEach, describe, expect, it, vi } from 'vitest';
import worker, { type Env } from '../src/index';

const files: Record<string, string> = {
  'training/latest/index.html': 'latest home',
  'training/latest/404.html': 'latest 404',
  'training/latest/hello_nextflow/index.html': 'hello',
  'training/latest/hello_nextflow.md': '# Hello',
  'training/latest/index.md': '# Home',
  'training/latest/pt.md': '# Início',
  'training/3.9/index.md': '# Old home',
  'training/latest/llms.txt': '# llms',
  'training/latest/_astro/app.abc123.js': 'js',
  'training/latest/pagefind/fragment/en_1.pf_fragment': 'bin',
  'training/development/index.html': 'dev home',
  'training/3.8.2/index.html': 'old home',
  'training/3.8.2/404.html': 'old 404',
  'training/3.8.2/es/a b/index.html': 'spaced',
  'training/versions.json': '[]',
};

const bucket = {
  async get(key: string) {
    return key in files
      ? { body: files[key], httpEtag: '"e"', httpMetadata: key.endsWith('.pf_fragment') ? { contentType: 'x/pf' } : {} }
      : null;
  },
  async head(key: string) {
    return key in files ? {} : null;
  },
};
const env = { BUCKET: bucket } as unknown as Env;

let cached: Map<string, Response>;
let pending: Promise<unknown>[];
const ctx = { waitUntil: (p: Promise<unknown>) => pending.push(p) } as unknown as ExecutionContext;

beforeEach(() => {
  cached = new Map();
  pending = [];
  vi.stubGlobal('caches', {
    default: {
      match: async (r: Request) => cached.get(r.url)?.clone(),
      put: async (r: Request, res: Response) => void cached.set(r.url, res),
    },
  });
});

const get = async (path: string, method = 'GET') => {
  const res = await worker.fetch(new Request(`https://training.nextflow.io${path}`, { method }), env, ctx);
  await Promise.all(pending);
  return res;
};

describe('router', () => {
  it('redirects / to /latest/, keeping the query', async () => {
    const res = await get('/?x=1');
    expect(res.status).toBe(302);
    expect(res.headers.get('Location')).toBe('https://training.nextflow.io/latest/?x=1');
  });

  it('serves index.html for directory paths', async () => {
    const res = await get('/latest/hello_nextflow/');
    expect(res.status).toBe(200);
    expect(await res.text()).toBe('hello');
    expect(res.headers.get('Content-Type')).toBe('text/html; charset=utf-8');
  });

  it('adds the trailing slash to directories, including dotted versions', async () => {
    expect((await get('/latest/hello_nextflow')).headers.get('Location')).toBe(
      'https://training.nextflow.io/latest/hello_nextflow/',
    );
    expect((await get('/3.8.2')).headers.get('Location')).toBe('https://training.nextflow.io/3.8.2/');
  });

  it("maps /<version>.md to the version's index.md", async () => {
    const latest = await get('/latest.md');
    expect(latest.status).toBe(200);
    expect(await latest.text()).toBe('# Home');
    expect(latest.headers.get('Content-Type')).toBe('text/markdown; charset=utf-8');
    expect(await (await get('/3.9.md')).text()).toBe('# Old home');
    const head = await get('/3.9.md', 'HEAD');
    expect(head.status).toBe(200);
    expect(head.body).toBeNull();
    expect(await (await get('/latest/pt.md')).text()).toBe('# Início');
    expect((await get('/4.0.md')).status).toBe(404);
  });

  it('decodes percent-encoded paths', async () => {
    expect(await (await get('/3.8.2/es/a%20b/')).text()).toBe('spaced');
  });

  it("returns the version's 404 page, falling back to latest", async () => {
    const old = await get('/3.8.2/nope/');
    expect(old.status).toBe(404);
    expect(await old.text()).toBe('old 404');
    const unknown = await get('/9.9/nope/');
    expect(unknown.status).toBe(404);
    expect(await unknown.text()).toBe('latest 404');
    expect(unknown.headers.get('Content-Type')).toBe('text/html; charset=utf-8');
    expect(cached.size).toBe(0);
  });

  it('applies the legacy redirects', async () => {
    const res = await get('/hello_nextflow/01_hello_world/?a=b');
    expect(res.status).toBe(301);
    expect(res.headers.get('Location')).toBe('https://training.nextflow.io/latest/hello_nextflow/01_hello_world/?a=b');
    expect((await get('/assets/x.png')).headers.get('Location')).toBe(
      'https://training.nextflow.io/latest/assets/x.png',
    );
  });

  it('sends noindex outside /latest/ only', async () => {
    expect((await get('/latest/')).headers.get('X-Robots-Tag')).toBeNull();
    expect((await get('/development/')).headers.get('X-Robots-Tag')).toBe('noindex');
    expect((await get('/3.8.2/')).headers.get('X-Robots-Tag')).toBe('noindex');
  });

  it('sets content types', async () => {
    const type = async (p: string) => (await get(p)).headers.get('Content-Type');
    expect(await type('/latest/hello_nextflow.md')).toBe('text/markdown; charset=utf-8');
    expect(await type('/latest/llms.txt')).toBe('text/plain; charset=utf-8');
    expect(await type('/versions.json')).toBe('application/json; charset=utf-8');
    expect(await type('/latest/pagefind/fragment/en_1.pf_fragment')).toBe('x/pf');
  });

  it('sets cache lifetimes and fills the cache', async () => {
    const cc = async (p: string) => (await get(p)).headers.get('Cache-Control');
    expect(await cc('/latest/_astro/app.abc123.js')).toBe('public, max-age=31536000, immutable');
    expect(await cc('/latest/')).toBe('public, max-age=300');
    expect(await cc('/development/')).toBe('public, max-age=300');
    expect(await cc('/versions.json')).toBe('public, max-age=300');
    expect(await cc('/3.8.2/')).toBe('public, max-age=86400');
    expect(cached.has('https://training.nextflow.io/3.8.2/')).toBe(true);
  });

  it('serves repeat requests from the cache', async () => {
    await get('/latest/');
    cached.set('https://training.nextflow.io/latest/', new Response('from cache'));
    expect(await (await get('/latest/')).text()).toBe('from cache');
  });

  it('answers HEAD like GET, without a body', async () => {
    const hit = await get('/3.8.2/', 'HEAD');
    expect(hit.status).toBe(200);
    expect(hit.body).toBeNull();
    expect(hit.headers.get('Content-Type')).toBe('text/html; charset=utf-8');
    expect(hit.headers.get('Cache-Control')).toBe('public, max-age=86400');
    const miss = await get('/3.8.2/nope/', 'HEAD');
    expect(miss.status).toBe(404);
    expect(miss.body).toBeNull();
    expect((await get('/3.8.2', 'HEAD')).headers.get('Location')).toBe('https://training.nextflow.io/3.8.2/');
    expect((await get('/', 'HEAD')).status).toBe(302);
    expect((await get('/hello_nextflow/', 'HEAD')).status).toBe(301);
    const versions = await get('/versions.json', 'HEAD');
    expect(versions.status).toBe(200);
    expect(versions.headers.get('Cache-Control')).toBe('public, max-age=300');
  });

  it('serves a cached GET to a later HEAD', async () => {
    await get('/latest/');
    const head = await get('/latest/', 'HEAD');
    expect(head.status).toBe(200);
    expect(head.body).toBeNull();
  });

  it('rejects other methods', async () => {
    expect((await get('/latest/', 'POST')).status).toBe(405);
  });
});
