#!/usr/bin/env node
/**
 * check-fences.mjs — static audit of Expressive Code meta strings.
 *
 * `astro build` does not validate meta. An off-by-one range, a label with
 * nowhere to sit, a marker pointing at code that scrolls out of view and a
 * plugin attribute with no plugin installed all build perfectly green. This
 * catches the mechanical subset of those.
 *
 *   node check-fences.mjs [dir]        # default: src/content/docs
 *   node check-fences.mjs --width 60   # override the marked-line budget
 *
 * Exit code 1 if any error-level finding is reported, so it can gate CI.
 * Warnings never fail the run.
 *
 * What it cannot check: whether the marked lines are the ones the prose is
 * talking about, and whether the rendered page looks right. Open a browser.
 */

import { readdirSync, readFileSync, statSync } from 'node:fs';
import { join, relative } from 'node:path';

const args = process.argv.slice(2);
const widthFlag = args.indexOf('--width');
const WIDTH = widthFlag === -1 ? 62 : Number(args[widthFlag + 1]);
const ROOT = args.find((a) => !a.startsWith('--') && a !== String(WIDTH)) ?? 'src/content/docs';

/** Attributes that silently do nothing unless their opt-in plugin is installed. */
const PLUGIN_ATTRS = [
  { re: /\bcollapse=\{/, pkg: '@expressive-code/plugin-collapsible-sections' },
  { re: /\bcollapseStyle=/, pkg: '@expressive-code/plugin-collapsible-sections' },
  { re: /\bshowLineNumbers\b/, pkg: '@expressive-code/plugin-line-numbers' },
  { re: /\bstartLineNumber=/, pkg: '@expressive-code/plugin-line-numbers' },
];

const findings = [];
const report = (level, file, line, message) => findings.push({ level, file, line, message });

/** Which opt-in plugins the project actually has, read from any package.json above us. */
function installedPlugins() {
  const found = new Set();
  for (const path of ['package.json', 'docs/package.json', '../package.json']) {
    try {
      const pkg = JSON.parse(readFileSync(path, 'utf8'));
      for (const dep of Object.keys({ ...pkg.dependencies, ...pkg.devDependencies })) found.add(dep);
    } catch {
      /* not there; fine */
    }
  }
  return found;
}

function* walk(dir) {
  for (const entry of readdirSync(dir)) {
    const path = join(dir, entry);
    if (statSync(path).isDirectory()) yield* walk(path);
    else if (/\.mdx?$/.test(path)) yield path;
  }
}

/** Parse `{1-3, 7}` / `ins={"Label:":5-8}` into ranges, ignoring inline "text" and /regex/ markers. */
function parseRanges(meta) {
  const ranges = [];
  for (const match of meta.matchAll(/(?:^|\s)(ins=|del=|mark=|collapse=)?\{([^}]*)\}/g)) {
    const [, kind = '', body] = match;
    const labelled = body.match(/^\s*"((?:[^"\\]|\\.)*)"\s*:\s*(.+)$/);
    const label = labelled?.[1];
    const spec = labelled ? labelled[2] : body;
    if (!/^[\d,\s-]+$/.test(spec)) continue;
    for (const part of spec.split(',')) {
      const trimmed = part.trim();
      if (!trimmed) continue;
      const [from, to] = trimmed.includes('-') ? trimmed.split('-').map(Number) : [Number(trimmed), Number(trimmed)];
      ranges.push({ kind: kind.replace('=', ''), from, to, label, collapse: kind === 'collapse=' });
    }
  }
  return ranges;
}

function checkFile(file, plugins) {
  const lines = readFileSync(file, 'utf8').split('\n');
  const where = (n) => `${relative(process.cwd(), file)}:${n}`;

  for (let i = 0; i < lines.length; i++) {
    const open = lines[i].match(/^(\s*)(`{3,})(\w[\w-]*)(.*)$/);
    if (!open) continue;
    const [, indent, ticks, lang, meta] = open;

    const body = [];
    let end = -1;
    for (let j = i + 1; j < lines.length; j++) {
      if (new RegExp(`^\\s{0,${indent.length}}${ticks}\\s*$`).test(lines[j])) {
        end = j;
        break;
      }
      body.push(lines[j]);
    }
    if (end === -1) continue; // unterminated; the build will complain
    i = end;
    if (!meta.trim()) continue;

    for (const { re, pkg } of PLUGIN_ATTRS) {
      if (re.test(meta) && !plugins.has(pkg)) {
        report('error', file, i + 1, `uses ${re.source.replace(/\\b|\\{/g, '')} but ${pkg} is not installed — silently ignored`);
      }
    }

    const ranges = parseRanges(meta);
    const collapsed = ranges.filter((r) => r.collapse);
    const marked = ranges.filter((r) => !r.collapse);

    for (const r of ranges) {
      if (r.from < 1 || r.to > body.length || r.from > r.to) {
        report('error', file, i + 1, `range {${r.from}-${r.to}} is outside the fence (${body.length} lines)`);
        continue;
      }
      // A long label needs a blank line of its own; otherwise it collides with code.
      if (r.label && r.label.length > 3 && body[r.from - 1].trim() !== '') {
        report('error', file, i + 1, `labelled marker "${r.label.slice(0, 32)}…" starts on line ${r.from}, which is not blank`);
      }
      if (r.collapse && r.to - r.from + 1 <= 5) {
        report('warn', file, i + 1, `collapse={${r.from}-${r.to}} hides only ${r.to - r.from + 1} lines — below the six-line floor`);
      }
    }

    for (const c of collapsed) {
      for (const m of marked) {
        if (m.from <= c.to && m.to >= c.from) {
          report('error', file, i + 1, `marker {${m.from}-${m.to}} is inside collapsed range {${c.from}-${c.to}}`);
        }
      }
    }

    // Width budget applies to marked, visible lines only.
    const hidden = (n) => collapsed.some((c) => n >= c.from && n <= c.to);
    const inlineMarkers = [...meta.matchAll(/(?:^|\s)(?:ins=|del=|mark=)?["']([^"']{2,})["']/g)].map((m) => m[1]);
    for (let n = 1; n <= body.length; n++) {
      if (hidden(n)) continue;
      const width = body[n - 1].length - indent.length;
      if (width <= WIDTH) continue;
      const isMarked =
        marked.some((r) => n >= r.from && n <= r.to) || inlineMarkers.some((text) => body[n - 1].includes(text));
      if (isMarked && !/\bwrap\b/.test(meta)) {
        report('warn', file, i + 1, `marked line ${n} is ${width} chars — the highlight may scroll out of view (add \`wrap\` or break the line)`);
      }
    }

    if (/\bins=|\bdel=/.test(meta) && lang !== 'diff') {
      report('warn', file, i + 1, `uses ins=/del= — confirm this block shows a real before/after, not a from-scratch example`);
    }
  }
}

const plugins = installedPlugins();
let files = 0;
try {
  for (const file of walk(ROOT)) {
    files++;
    checkFile(file, plugins);
  }
} catch (error) {
  console.error(`check-fences: cannot read ${ROOT} — ${error.message}`);
  process.exit(2);
}

const errors = findings.filter((f) => f.level === 'error');
const warnings = findings.filter((f) => f.level === 'warn');

for (const f of [...errors, ...warnings]) {
  const tag = f.level === 'error' ? 'ERROR' : 'warn ';
  console.log(`${tag} ${relative(process.cwd(), f.file)}:${f.line}  ${f.message}`);
}

console.log(
  `\nchecked ${files} file(s): ${errors.length} error(s), ${warnings.length} warning(s)` +
    (findings.length ? '' : ' — all meta strings consistent'),
);
console.log('This checks the mechanical parts only. Whether the right lines are marked, and whether the');
console.log('page reads well, needs a browser.');

process.exit(errors.length ? 1 : 0);
