/**
 * Sätteri hast plugin: a "Takeaway" heading and everything after it, up to the next heading or `---`,
 * become one `<section class="nf-takeaway">` (styles/content.css).
 *
 * English pages match the heading text. A translated page boxes the headings at the same positions as
 * the English page's "Takeaway" headings, so it needs no translated words and no heading levels. A
 * translation with a different number of headings (out of date) gets no boxes until it is re-translated.
 */
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { translatedLocales } from '../locales.mjs';

const TITLE = 'Takeaway';
const HEADING = /^h[1-6]$/;
const BOUNDARY = /^(h[1-6]|hr)$/;
const isElement = (node, test) => node?.type === 'element' && test.test(node.tagName);

/** Heading texts of a Markdown/MDX source, in order, skipping frontmatter and code. */
export function sourceHeadings(text) {
  const lines = text.split('\n');
  let i = lines[0]?.trim() === '---' ? lines.indexOf('---', 1) + 1 : 0;
  const out = [];
  let fence = '';
  for (; i < lines.length; i++) {
    const line = lines[i];
    const f = /^\s*(`{3,}|~{3,})/.exec(line);
    if (fence) {
      if (f && f[1][0] === fence[0] && f[1].length >= fence.length) fence = '';
      continue;
    }
    if (f) fence = f[1];
    const h = !f && /^ {0,3}#{1,6}\s+(.*?)\s*(?:\{#[^}]*\})?\s*#*\s*$/.exec(line);
    if (h) out.push(h[1]);
  }
  return out;
}

const englishFlags = new Map();
/** For a translated page, which heading positions are takeaways in English. Undefined for English. */
function translatedFlags(fileURL) {
  if (!fileURL) return;
  const path = fileURLToPath(fileURL);
  const m = new RegExp(`([/\\\\]content[/\\\\]docs[/\\\\])(${translatedLocales.join('|')})[/\\\\]`).exec(path);
  if (!m) return;
  const english = path.slice(0, m.index) + m[1] + path.slice(m.index + m[0].length);
  if (!englishFlags.has(english)) {
    let flags = null;
    try {
      flags = sourceHeadings(readFileSync(english, 'utf8')).map((t) => t === TITLE);
    } catch {}
    englishFlags.set(english, flags);
  }
  return englishFlags.get(english);
}

export function takeaways() {
  return {
    name: 'nf-takeaways',
    after(root, ctx) {
      // Every heading in document order, with the top-level child it sits in.
      const headings = [];
      const walk = (node, top) => {
        if (isElement(node, HEADING)) headings.push({ node, top });
        for (const child of node.children ?? []) walk(child, top);
      };
      root.children.forEach((child, i) => walk(child, i));

      const flags = translatedFlags(ctx.fileURL);
      const isTakeaway =
        flags === undefined
          ? (h) => ctx.textContent(h.node).trim() === TITLE
          : flags && flags.length === headings.length
            ? (h, i) => flags[i]
            : () => false;
      const starts = new Set(
        headings.filter((h, i) => isTakeaway(h, i) && root.children[h.top] === h.node).map((h) => h.top),
      );
      if (!starts.size) return;

      const children = root.children;
      const out = [];
      for (let i = 0; i < children.length; i++) {
        if (!starts.has(i)) {
          out.push(children[i]);
          continue;
        }
        const box = [children[i]];
        while (i + 1 < children.length && !isElement(children[i + 1], BOUNDARY)) box.push(children[++i]);
        out.push({ type: 'element', tagName: 'section', properties: { className: ['nf-takeaway'] }, children: box });
      }
      ctx.setField(root, 'children', out);
    },
  };
}
