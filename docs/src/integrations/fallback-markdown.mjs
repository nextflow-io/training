import { copyFile, mkdir, readdir } from 'node:fs/promises';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

async function* markdownFiles(dir, root = dir) {
  for (const entry of await readdir(dir, { withFileTypes: true })) {
    const p = path.join(dir, entry.name);
    if (entry.isDirectory()) yield* markdownFiles(p, root);
    else if (entry.name.endsWith('.md')) yield path.relative(root, p);
  }
}

/**
 * starlight-page-actions writes a raw .md only for pages that have a source file, but
 * fallback pages (English content under a locale) link to one too: give them the English copy.
 */
export default function fallbackMarkdown(locales) {
  return {
    name: 'fallback-markdown',
    hooks: {
      'astro:build:done': async ({ dir, logger }) => {
        const dist = fileURLToPath(dir);
        const english = [];
        for await (const file of markdownFiles(dist)) {
          if (!locales.includes(file.split(path.sep)[0].replace(/\.md$/, ''))) english.push(file);
        }
        let copied = 0;
        for (const locale of locales) {
          for (const file of english) {
            // The page actions ask for `/<base>/<locale>.md` on a locale's home page.
            const target = file === 'index.md' ? path.join(dist, `${locale}.md`) : path.join(dist, locale, file);
            try {
              await copyFile(path.join(dist, file), target, 1); // COPYFILE_EXCL: never overwrite a translation
              copied++;
            } catch (e) {
              if (e.code === 'ENOENT') {
                await mkdir(path.dirname(target), { recursive: true });
                await copyFile(path.join(dist, file), target, 1);
                copied++;
              } else if (e.code !== 'EEXIST') throw e;
            }
          }
        }
        logger.info(`${copied} fallback Markdown files`);
      },
    },
  };
}
