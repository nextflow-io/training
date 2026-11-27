// Small helpers shared by the client-side views of My courses.
import { icon } from '../icons.mjs';

export const esc = (s: string) =>
  String(s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]!);

/** UI strings a component embedded as JSON in `data-t`, with `{{name}}` placeholders. */
export function strings(el: HTMLElement) {
  const all: Record<string, string> = JSON.parse(el.dataset.t || '{}');
  return (key: string, vars: Record<string, string | number> = {}) =>
    (all[key] ?? key).replace(/\{\{(\w+)\}\}/g, (_, k) => String(vars[k] ?? ''));
}

/** A progress ring. `f` runs from 0 to 1. */
export function ring(f: number, size: number, width: number) {
  const r = (size - width) / 2;
  const c = 2 * Math.PI * r;
  const mid = size / 2;
  const arc =
    f > 0
      ? `<circle cx="${mid}" cy="${mid}" r="${r}" fill="none" stroke="var(--nf-green)" stroke-width="${width}" stroke-linecap="round" stroke-dasharray="${c.toFixed(2)}" stroke-dashoffset="${(c * (1 - Math.min(f, 1))).toFixed(2)}" transform="rotate(-90 ${mid} ${mid})"/>`
      : '';
  return `<svg class="nf-ring" width="${size}" height="${size}" viewBox="0 0 ${size} ${size}" aria-hidden="true"><circle cx="${mid}" cy="${mid}" r="${r}" fill="none" stroke="var(--nf-line)" stroke-width="${width}"/>${arc}</svg>`;
}

/** The 40px ring with the percentage inside, or a filled tick when complete. */
export function courseRing(f: number) {
  const inner =
    f >= 1 ? icon('check', 14, 3, 'nf-ring-mark') : `<span class="nf-ring-pct">${Math.round(f * 100)}</span>`;
  return `<span class="nf-course-ring${f >= 1 ? ' is-done' : ''}">${ring(f, 40, 4)}${inner}</span>`;
}

/** Announce a change to screen readers through the header's live region. */
export function announce(msg: string) {
  const el = document.getElementById('nf-announce');
  if (!el) return;
  el.textContent = '';
  requestAnimationFrame(() => (el.textContent = msg));
}

export const reducedMotion = () => matchMedia('(prefers-reduced-motion: reduce)').matches;
