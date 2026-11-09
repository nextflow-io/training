import assert from "node:assert/strict";
import { mkdirSync, mkdtempSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { test } from "node:test";
import { changedPagesMarkdown } from "./preview-pages.mjs";

const URL = "https://pr-7.example.workers.dev";
const DOCS = "docs/src/content/docs/";

// An English-only build, plus one German page, a redirect left at a deleted page's path and
// pages with hostile names, since a fork controls the artifact too.
const HOSTILE = ["a](javascript:alert(1))[b", "a`b", "<img src=x>", "héllo"];
const site = mkdtempSync(join(tmpdir(), "preview-pages-"));
const pages = [
  "",
  "hello_nextflow",
  "hello_nextflow/01_hello_world",
  "old_page",
  "de/hello_nextflow/01_hello_world",
  ...HOSTILE,
];
for (let i = 1; i <= 60; i++) pages.push(`many/p${String(i).padStart(2, "0")}`);
for (const page of pages) {
  mkdirSync(join(site, page), { recursive: true });
  writeFileSync(join(site, page, "index.html"), "");
}

const changed = (...names) =>
  names.map((filename) => ({ filename, status: "modified" }));
const run = (files) => changedPagesMarkdown(files, site, URL);

test("links an edited English page", () => {
  assert.equal(
    run(changed(`${DOCS}hello_nextflow/01_hello_world.mdx`)),
    `**Changed pages**\n\n- [\`/hello_nextflow/01_hello_world/\`](${URL}/latest/hello_nextflow/01_hello_world/)`,
  );
});

test("maps index.mdx to its folder", () => {
  const md = run(
    changed(`${DOCS}index.mdx`, `${DOCS}hello_nextflow/index.mdx`),
  );
  assert.match(
    md,
    /^- \[`\/`\]\(https:\/\/pr-7\.example\.workers\.dev\/latest\/\)$/m,
  );
  assert.match(
    md,
    /^- \[`\/hello_nextflow\/`\]\(.*\/latest\/hello_nextflow\/\)$/m,
  );
});

test("skips pages missing from the build, such as translations in an English-only build", () => {
  assert.equal(
    run(
      changed(`${DOCS}fr/hello_nextflow/01_hello_world.mdx`, `${DOCS}nope.mdx`),
    ),
    "",
  );
});

test("groups translations under the page, English first", () => {
  const md = run(
    changed(`${DOCS}de/hello_nextflow/01_hello_world.mdx`, `${DOCS}index.mdx`),
  );
  assert.deepEqual(md.split("\n").slice(2), [
    `- [\`/\`](${URL}/latest/)`,
    `- \`/hello_nextflow/01_hello_world/\` · [de](${URL}/latest/de/hello_nextflow/01_hello_world/)`,
  ]);
  assert.match(
    run(
      changed(
        `${DOCS}hello_nextflow/01_hello_world.mdx`,
        `${DOCS}de/hello_nextflow/01_hello_world.mdx`,
      ),
    ),
    /^- \[`\/hello_nextflow\/01_hello_world\/`\]\(\S+\) · \[de\]\(\S+\)$/m,
  );
});

test("uses the new path of a rename and skips deletions", () => {
  const md = run([
    {
      filename: `${DOCS}hello_nextflow/01_hello_world.mdx`,
      previous_filename: `${DOCS}x.mdx`,
      status: "renamed",
    },
    { filename: `${DOCS}old_page.mdx`, status: "removed" },
  ]);
  assert.match(md, /hello_nextflow\/01_hello_world/);
  assert.doesNotMatch(md, /old_page/);
});

test("caps the list", () => {
  const many = (n) =>
    changed(
      ...Array.from(
        { length: n },
        (_, i) => `${DOCS}many/p${String(i + 1).padStart(2, "0")}.mdx`,
      ),
    );
  const rows = (md) => md.split("\n").filter((l) => l.startsWith("- ")).length;
  assert.equal(rows(run(many(1))), 1);
  assert.equal(rows(run(many(5))), 5);
  const md = run(many(60));
  assert.equal(rows(md), 10);
  assert.match(md, /\/many\/p10\//);
  assert.match(md, /\n\n…and 50 more\.$/);
});

test("notes when the API truncated the file list", () => {
  const files = [
    ...changed(`${DOCS}index.mdx`),
    ...changed(...Array(2999).fill("nextflow-run/x.nf")),
  ];
  assert.match(run(files), /first 3,000 changed files/);
  assert.doesNotMatch(run(files.slice(0, -1)), /3,000/);
});

test("ignores hostile and non-page file names", () => {
  for (const name of [
    ...HOSTILE.map((page) => `${DOCS}${page}.mdx`),
    `${DOCS}../../../hello_nextflow.mdx`,
    `${DOCS}hello_nextflow/01_hello_world.mdx\n- [x](https://evil.example)`,
    `${DOCS}hello_nextflow/01_hello_world.md`,
    "x/docs/src/content/docs/index.mdx",
    "hello-nextflow/hello-world.nf",
  ]) {
    assert.equal(run(changed(name)), "", name);
  }
});

test("says nothing when no page changed", () => {
  assert.equal(run([]), "");
  assert.equal(run(changed(".github/workflows/docs-preview.yml")), "");
});
