import { existsSync } from "node:fs";
import { join } from "node:path";

// File names come from the PR, so only these characters may reach the comment.
// Two-letter top folders are translations, as in docs-build.yml.
const PAGE =
  /^docs\/src\/content\/docs\/(?:([a-z]{2})\/)?((?:[\w-]+\/)*[\w-]+)\.mdx$/;
const MAX_ROWS = 10;

/**
 * Markdown that links the preview of each changed page, or '' if there are none.
 * `files` is the pulls.listFiles response; a page counts only if `siteDir` has it.
 */
export function changedPagesMarkdown(files, siteDir, previewUrl) {
  const pages = new Map();
  for (const { filename, status } of files) {
    const match = status !== "removed" && PAGE.exec(filename);
    if (!match) continue;
    const [, lang = "", path] = match;
    const slug = path.toLowerCase().replace(/(^|\/)index$/, "");
    if (!existsSync(join(siteDir, lang, slug, "index.html"))) continue;
    if (!pages.has(slug)) pages.set(slug, new Set());
    pages.get(slug).add(lang);
  }
  if (!pages.size) return "";

  const base = `${previewUrl.replace(/\/+$/, "")}/latest/`;
  const href = (lang, slug) =>
    base + (lang && `${lang}/`) + (slug && `${slug}/`);
  const rows = [...pages]
    .sort(([a, x], [b, y]) => y.has("") - x.has("") || a.localeCompare(b))
    .map(([slug, langs]) => {
      const name = `\`/${slug}${slug ? "/" : ""}\``;
      const translations = [...langs]
        .filter(Boolean)
        .sort()
        .map((lang) => `[${lang}](${href(lang, slug)})`);
      const english = langs.has("") ? `[${name}](${href("", slug)})` : name;
      return `- ${[english, translations.join(" ")].filter(Boolean).join(" · ")}`;
    });

  const lines = ["**Changed pages**", "", ...rows.slice(0, MAX_ROWS)];
  if (rows.length > MAX_ROWS)
    lines.push("", `…and ${rows.length - MAX_ROWS} more.`);
  if (files.length >= 3000) {
    lines.push(
      "",
      "GitHub lists only the first 3,000 changed files, so some pages may be missing.",
    );
  }
  return lines.join("\n");
}
