import path from 'node:path';
import { fileURLToPath } from 'node:url';

const CONTENT_DIR = '/src/content/docs/';
const SCHEME = /^[a-z][a-z\d+.-]*:/i;
const MD_LINK = /^(?<file>[^?#]*?)\.mdx?(?<rest>[?#].*)?$/;
const LINK_COMPONENTS = ['LinkButton', 'LinkCard'];

/** Rewrites one href from a content file to a site URL, or returns undefined to leave it as is. */
export function resolveHref(href, file, base) {
  if (typeof href !== 'string' || !href || href.startsWith('#') || href.startsWith('//') || SCHEME.test(href)) return;
  const root = base.endsWith('/') ? base : `${base}/`;
  if (href.startsWith('/')) return root === '/' || href.startsWith(root) ? undefined : root + href.slice(1);
  const m = MD_LINK.exec(href);
  if (!m) return;
  const target = path.posix.normalize(path.posix.join(path.posix.dirname(file), m.groups.file));
  // Starlight slugs are lowercase, and `index` is the directory URL.
  const slug = target.toLowerCase().replace(/(^|\/)index$/, '');
  return root + (slug && slug !== '.' ? `${slug}/` : '') + (m.groups.rest ?? '');
}

function contentPath(fileURL) {
  if (!fileURL) return;
  const p = fileURLToPath(fileURL).replaceAll('\\', '/');
  const i = p.lastIndexOf(CONTENT_DIR);
  return i === -1 ? undefined : p.slice(i + CONTENT_DIR.length);
}

function jsxHref(node, ctx, base) {
  const file = contentPath(ctx.fileURL);
  if (!file) return;
  let changed = false;
  const attributes = node.attributes.map((a) => {
    if (a.type !== 'mdxJsxAttribute' || a.name !== 'href') return a;
    const href = resolveHref(a.value, file, base);
    if (href === undefined) return a;
    changed = true;
    return { ...a, value: href };
  });
  // setProperty does not reach JSX attributes; a returned node replaces this one.
  if (changed) return { ...node, attributes };
}

/** Sätteri hast plugin: content links resolve to `<base><locale>/<slug>/`, and root-absolute links get the base. */
export default function links({ base = '/' } = {}) {
  return {
    name: 'training-links',
    element: {
      filter: ['a'],
      visit(node, ctx) {
        const file = contentPath(ctx.fileURL);
        const href = file && resolveHref(node.properties?.href, file, base);
        if (href !== undefined && file) ctx.setProperty(node, 'href', href);
      },
    },
    mdxJsxFlowElement: { filter: LINK_COMPONENTS, visit: (node, ctx) => jsxHref(node, ctx, base) },
    mdxJsxTextElement: { filter: LINK_COMPONENTS, visit: (node, ctx) => jsxHref(node, ctx, base) },
  };
}
