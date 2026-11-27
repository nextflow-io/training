/**
 * Sätteri plugin: `:::type[Title]` admonitions in every MkDocs Material type, optionally collapsible.
 *
 *   :::exercise[Title]            a static box
 *   :::solution[Title]{collapsed} collapsible, closed (MkDocs `???`)
 *   :::tip[Title]{open}           collapsible, open (MkDocs `???+`)
 *
 * Static note, tip, caution and danger are left to Starlight, which translates their default titles.
 * Everything else renders here. Colours and icons are in styles/custom.css.
 */

export const types = [
  'note',
  'abstract',
  'info',
  'tip',
  'success',
  'question',
  'warning',
  'caution',
  'failure',
  'danger',
  'bug',
  'example',
  'quote',
  'exercise',
  'solution',
  'full-code',
  'optional',
  'courses',
  'catalog',
];
const known = new Set(types);
const starlight = new Set(['note', 'tip', 'caution', 'danger']);

/** A default title: `full-code` → "Full code". */
export const defaultTitle = (type) => (type.charAt(0).toUpperCase() + type.slice(1)).replace(/-/g, ' ');

const element = (tagName, properties, children = []) => ({
  type: 'paragraph',
  data: { hName: tagName, hProperties: properties },
  children,
});

export function admonitions() {
  return {
    name: 'nf-admonitions',
    containerDirective(node, ctx) {
      if (!known.has(node.name)) return;
      const attrs = node.attributes ?? {};
      const collapsible = 'collapsed' in attrs || 'open' in attrs;
      if (!collapsible && starlight.has(node.name)) return;

      const children = [...node.children];
      let titleNodes = [{ type: 'text', value: defaultTitle(node.name) }];
      let title = defaultTitle(node.name);
      const first = children[0];
      if (first?.type === 'paragraph' && first.data?.directiveLabel && first.children.length > 0) {
        titleNodes = first.children;
        title = ctx.textContent(first);
        children.shift();
      }
      const className = ['nf-admonition', `nf-admonition--${node.name}`, ...(attrs.class ?? '').split(/\s+/)]
        .filter(Boolean)
        .join(' ');
      const content = element('div', { class: 'nf-admonition__content' }, children);

      if (collapsible) {
        return element('details', { class: className, ...('open' in attrs && { open: true }) }, [
          element('summary', { class: 'nf-admonition__title' }, titleNodes),
          content,
        ]);
      }
      return element('aside', { class: className, 'aria-label': title }, [
        element('p', { class: 'nf-admonition__title', 'aria-hidden': 'true' }, titleNodes),
        content,
      ]);
    },
  };
}
