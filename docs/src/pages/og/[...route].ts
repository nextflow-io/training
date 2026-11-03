import { getCollection } from 'astro:content';
import { OGImageRoute } from 'astro-og-canvas';
import config from 'virtual:starlight/user-config';

const entries = await getCollection('docs');
const titles = new Map(entries.map((entry) => [entry.id, entry.data.title]));
const homes = new Set(['index', ...Object.keys(config.locales)]);

const families = ['Inter', 'Noto Sans KR', 'Noto Sans Devanagari'];
const green: [number, number, number] = [13, 192, 157];
// Word joiners stop Korean wrapping mid-word, like CSS `word-break: keep-all`.
const keepAll = (text: string) => text.replace(/(?<=\p{Script=Hangul})(?=\p{Script=Hangul})/gu, '\u2060');

export const { getStaticPaths, GET } = await OGImageRoute({
  pages: Object.fromEntries(entries.map((entry) => [entry.id, entry.data])),
  // The default slug drops everything after the last dot in the id.
  getSlug: (id) => `${id}.png`,
  getImageOptions: (id, page) => {
    const parent = id.split('/').slice(0, -1).join('/');
    return {
      // Home pages are titled "Home"; their description says more.
      title: keepAll(homes.has(id) ? (page.description ?? page.title) : page.title),
      description: homes.has(parent) ? '' : (titles.get(parent) ?? ''),
      logo: { path: 'src/assets/nextflow-logo-white.png', size: [400] },
      // A gradient makes each card about 10 ms slower to encode.
      bgGradient: [[22, 15, 38]],
      border: { color: green, width: 16, side: 'inline-start' },
      padding: 64,
      fonts: [
        'src/assets/og/Inter-Bold.ttf',
        'src/assets/og/NotoSansKR-Bold.ttf',
        'src/assets/og/NotoSansDevanagari-Bold.ttf',
      ],
      font: {
        title: { families, weight: 'Bold', size: 64, lineHeight: 1.2, color: [255, 255, 255] },
        description: { families, weight: 'Bold', size: 36, color: green },
      },
    };
  },
});
