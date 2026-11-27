import assert from 'node:assert/strict';
import { test } from 'node:test';
import { createStore, formatMinutes, merge, validate } from './my-courses.ts';

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

test('logging out clears the browser copy', async () => {
  const backend = memory();
  const store = createStore(items, backend);
  store.add('nextflow_run');
  store.setParts('hello_nextflow', ['hello_nextflow/01_a'], true);
  await store.logout();
  assert.equal(backend.saved, null);
  assert.deepEqual(store.basket(), []);
  assert.equal(store.done('hello_nextflow'), 0);
});

test('an older version keeps courses and parts it does not know', () => {
  const backend = memory({
    basket: ['hello_nextflow', 'new_course'],
    progress: { new_course: ['new_course/01'], hello_nextflow: ['hello_nextflow/01_a', 'hello_nextflow/09_new'] },
  });
  const store = createStore(items, backend);
  assert.deepEqual(store.basket(), ['hello_nextflow']);
  store.add('nextflow_run');
  store.setParts('hello_nextflow', ['hello_nextflow/01_a'], false);
  store.remove('hello_nextflow');
  assert.deepEqual(backend.saved.basket, ['nextflow_run', 'new_course']);
  assert.deepEqual(backend.saved.progress, {
    new_course: ['new_course/01'],
    hello_nextflow: ['hello_nextflow/09_new'],
  });
});

test('merge is a union and keeps the account goal', () => {
  const s = merge(
    { goal: 'write', basket: ['a'], progress: { a: ['a/1'] } },
    { goal: 'custom', basket: ['b', 'a'], progress: { a: ['a/2', 'a/1'], b: ['b/1'] } },
  );
  assert.deepEqual(s, { version: 1, goal: 'write', basket: ['a', 'b'], progress: { a: ['a/1', 'a/2'], b: ['b/1'] } });
  assert.equal(merge(null, { goal: 'custom' }).goal, 'custom');
});

const ada = { sub: 'seqera|ada', name: 'Ada Lovelace', email: 'ada@example.com' };
const remote = (server) => ({
  server,
  saves: 0,
  loggedOut: false,
  async load() {
    return this.loggedOut ? null : { user: ada, state: this.server };
  },
  async save(s) {
    if (this.loggedOut) return false;
    this.saves++;
    this.server = JSON.parse(JSON.stringify(s));
    return true;
  },
  logout() {},
});

test('the browser copy merges into the account once', async () => {
  const backend = memory({ basket: ['nextflow_run'], progress: { nextflow_run: ['nextflow_run/01_a'] } });
  const account = remote({
    goal: 'write',
    basket: ['hello_nextflow'],
    progress: { hello_nextflow: ['hello_nextflow/01_a'] },
  });
  const store = createStore(items, backend, account);
  assert.equal(store.user, null);
  await store.ready;
  assert.deepEqual(store.user, { name: 'Ada Lovelace', email: 'ada@example.com', initials: 'AL' });
  assert.deepEqual(store.basket(), ['nextflow_run', 'hello_nextflow']);
  assert.equal(store.state.goal, 'write');
  assert.deepEqual(account.server.basket, ['nextflow_run', 'hello_nextflow']);
  assert.equal(backend.saved.account.sub, 'seqera|ada');
  assert.equal(backend.saved.pending, undefined);

  // Another device un-marks a part: the next page load takes the account as it is.
  account.server = { ...account.server, progress: { hello_nextflow: ['hello_nextflow/01_a'] } };
  const next = createStore(items, backend, account);
  assert.deepEqual(next.user, store.user);
  await next.ready;
  assert.equal(next.done('nextflow_run'), 0);
  assert.equal(account.saves, 1);
});

test('changes are pushed, and unsaved ones merge on the next load', async () => {
  const account = remote(null);
  const backend = memory();
  const store = createStore(items, backend, account);
  store.add('nextflow_run');
  await store.ready;
  assert.deepEqual(account.server.basket, ['nextflow_run']);
  store.setParts('hello_nextflow', ['hello_nextflow/01_a'], true);
  assert.equal(backend.saved.pending, true);
  await store.flush();
  assert.deepEqual(account.server.progress, { hello_nextflow: ['hello_nextflow/01_a'] });
  assert.equal(backend.saved.pending, undefined);

  // Offline when the change was made: the next load merges it in rather than dropping it.
  const save = account.save;
  account.save = async () => {
    throw new Error('offline');
  };
  store.setParts('nextflow_run', ['nextflow_run/02_b'], true);
  await store.flush();
  account.save = save;
  const next = createStore(items, backend, account);
  await next.ready;
  assert.deepEqual(account.server.progress.nextflow_run, ['nextflow_run/02_b']);
});

test("another account's copy is replaced, not merged", async () => {
  const backend = memory({ account: { ...ada, sub: 'seqera|bob' }, basket: ['nextflow_run'], progress: {} });
  const account = remote({ basket: [], progress: {} });
  const store = createStore(items, backend, account);
  await store.ready;
  assert.deepEqual(store.basket(), []);
  assert.equal(account.saves, 0);
});

test('an expired session drops the account copy and keeps anonymous progress', async () => {
  const tagged = memory({ account: ada, basket: ['nextflow_run'], progress: {} });
  const expired = remote(null);
  expired.loggedOut = true;
  const a = createStore(items, tagged, expired);
  await a.ready;
  assert.equal(a.user, null);
  assert.equal(tagged.saved, null);

  const anonymous = memory({ basket: ['nextflow_run'], progress: {} });
  const b = createStore(items, anonymous, expired);
  await b.ready;
  assert.deepEqual(b.basket(), ['nextflow_run']);
});

test('formatMinutes', () => {
  assert.equal(formatMinutes(65), '1 h 5 min');
  assert.equal(formatMinutes(120), '2 h');
  assert.equal(formatMinutes(45), '45 min');
});
