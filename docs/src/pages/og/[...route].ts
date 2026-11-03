import { getCollection } from 'astro:content';
import { OGImageRoute } from 'astro-og-canvas';

const entries = await getCollection('docs');

export const { getStaticPaths, GET } = await OGImageRoute({
  pages: Object.fromEntries(entries.map((entry) => [entry.id, entry.data])),
  getImageOptions: (_path, page) => ({
    title: page.title,
    description: page.description ?? '',
    bgGradient: [[11, 25, 34]],
    border: { color: [13, 192, 157], width: 16, side: 'inline-start' },
    padding: 60,
  }),
});
