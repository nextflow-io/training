import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createStore, formatMinutes, validate } from './my-courses.ts';

const items = [
  { id: 'nextflow_run', title: 'Nextflow Run', parts: ['nextflow_run/01_a', 'nextflow_run/02_b'], minutes: 45 },
  { id: 'hello_nextflow', title: 'Hello Nextflow', parts: ['hello_nextflow/01_a'], minutes: 30 },
  { id: 'side_quests/metadata', title: 'Metadata', parts: ['side_quests/metadata'], minutes: 45, sq: true },
];
const memory = (initial = null) => ({
  saved: initial,
  load() {
    return this.saved;
  },
  save(s) {
    this.saved = s;
  },
  clear() {
    this.saved = null;
  },
});

test('validate drops unknown ids and parts and orders the basket', () => {
  const s = validate(
    {
      goal: 'nope',
      basket: ['side_quests/metadata', 'gone', 'nextflow_run'],
      progress: { nextflow_run: ['x', 'nextflow_run/02_b'], gone: ['a'] },
    },
    items,
  );
  assert.deepEqual(s, {
    version: 1,
    goal: 'run',
    basket: ['nextflow_run', 'side_quests/metadata'],
    progress: { nextflow_run: ['nextflow_run/02_b'] },
  });
  assert.deepEqual(validate('garbage', items).basket, []);
});

test('marking a part done adds the item; removing keeps progress', () => {
  const backend = memory();
  const store = createStore(items, backend);
  const changes = [];
  store.subscribe((c) => changes.push(c.type));
  store.setParts('hello_nextflow', ['hello_nextflow/01_a'], true);
  assert.ok(store.has('hello_nextflow'));
  assert.equal(store.fraction('hello_nextflow'), 1);
  assert.deepEqual(changes, ['add', 'progress']);
  store.add('nextflow_run');
  assert.deepEqual(store.basket(), ['nextflow_run', 'hello_nextflow']);
  assert.equal(store.overall(), 1 / 3);
  store.remove('hello_nextflow');
  assert.deepEqual(store.hidden(), ['hello_nextflow']);
  store.reset();
  assert.deepEqual(store.basket(), []);
  assert.equal(store.done('hello_nextflow'), 1);
  store.setParts('hello_nextflow', ['hello_nextflow/01_a'], false);
  assert.deepEqual(backend.saved.progress, {});
});

test('logging out clears the browser copy', () => {
  const backend = memory();
  const store = createStore(items, backend);
  store.add('nextflow_run');
  store.setParts('hello_nextflow', ['hello_nextflow/01_a'], true);
  store.logout();
  assert.equal(backend.saved, null);
  assert.deepEqual(store.basket(), []);
  assert.equal(store.done('hello_nextflow'), 0);
});

test('formatMinutes', () => {
  assert.equal(formatMinutes(65), '1 h 5 min');
  assert.equal(formatMinutes(120), '2 h');
  assert.equal(formatMinutes(45), '45 min');
});
