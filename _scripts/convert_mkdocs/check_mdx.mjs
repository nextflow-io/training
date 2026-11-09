// Compile every .mdx file under a directory with Sätteri, the parser `astro build` uses,
// and report syntax errors. Seconds for the whole corpus, instead of a full build.
//
//   node _scripts/convert_mkdocs/check_mdx.mjs <dir> [--quiet]
//
// Needs `npm ci` in docs/ first (it loads satteri from docs/node_modules).
import { readFileSync, readdirSync, statSync } from "node:fs";
import { createRequire } from "node:module";
import path from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const require = createRequire(path.join(here, "../../docs/package.json"));
const { mdxToJs } = await import(
  pathToFileURL(require.resolve("satteri")).href
);

const [dir, ...flags] = process.argv.slice(2);
if (!dir) {
  console.error("usage: check_mdx.mjs <dir> [--quiet]");
  process.exit(2);
}

function* walk(d) {
  for (const name of readdirSync(d).sort()) {
    const p = path.join(d, name);
    if (statSync(p).isDirectory()) yield* walk(p);
    else if (p.endsWith(".mdx")) yield p;
  }
}

let files = 0;
let failed = 0;
for (const file of walk(dir)) {
  files++;
  const src = readFileSync(file, "utf8").replace(
    /^---\n[\s\S]*?\n---\n/,
    (fm) => "\n".repeat(fm.split("\n").length - 1),
  );
  try {
    mdxToJs(src, {
      features: { directive: true, gfm: true, headingAttributes: true },
    });
  } catch (e) {
    failed++;
    console.log(
      `${path.relative(dir, file)}: ${String(e.message ?? e).split("\n")[0]}`,
    );
  }
}
if (!flags.includes("--quiet") || failed)
  console.log(`${files} files, ${failed} with MDX errors`);
process.exit(failed ? 1 : 0);
