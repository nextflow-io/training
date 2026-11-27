import assert from 'node:assert/strict';
import { test } from 'node:test';
import { sourceHeadings } from './takeaways.mjs';

test('sourceHeadings skips frontmatter and code', () => {
  const src = [
    '---',
    'title: Page',
    '---',
    '## 1. Section {#one}',
    '```bash',
    '# not a heading',
    '```',
    '### Takeaway',
    ':::tip',
    '#### Inside an aside',
    ':::',
    '### What’s next?',
  ].join('\n');
  assert.deepEqual(sourceHeadings(src), ['1. Section', 'Takeaway', 'Inside an aside', 'What’s next?']);
});
