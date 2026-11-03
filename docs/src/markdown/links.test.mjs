import assert from 'node:assert/strict';
import { test } from 'node:test';
import { resolveHref } from './links.mjs';

const page = 'hello_nextflow/01_hello_world.mdx';
const ptPage = 'pt/hello_nextflow/01_hello_world.mdx';

test('relative .md and .mdx links resolve to the page URL under the base', () => {
  assert.equal(resolveHref('00_orientation.md', page, '/latest/'), '/latest/hello_nextflow/00_orientation/');
  assert.equal(
    resolveHref('./02_hello_channels.mdx#21-foo', page, '/latest/'),
    '/latest/hello_nextflow/02_hello_channels/#21-foo',
  );
  assert.equal(resolveHref('../envsetup/index.md', page, '/3.9/'), '/3.9/envsetup/');
  assert.equal(resolveHref('../index.md', page, '/'), '/');
  assert.equal(resolveHref('README.md', page, '/'), '/hello_nextflow/readme/');
});

test('the locale comes from the page path', () => {
  assert.equal(resolveHref('../envsetup/index.md', ptPage, '/latest/'), '/latest/pt/envsetup/');
  assert.equal(resolveHref('../index.md', ptPage, '/latest/'), '/latest/pt/');
});

test('root-absolute links get the base once', () => {
  assert.equal(resolveHref('/envsetup/', page, '/latest/'), '/latest/envsetup/');
  assert.equal(resolveHref('/latest/envsetup/', page, '/latest/'), undefined);
  assert.equal(resolveHref('/envsetup/', page, '/'), undefined);
});

test('other links are left alone', () => {
  for (const href of [
    '#anchor',
    'https://nextflow.io/docs/',
    'mailto:a@b.c',
    '//cdn.example/x',
    'img/x.png',
    undefined,
  ]) {
    assert.equal(resolveHref(href, page, '/latest/'), undefined, String(href));
  }
});
