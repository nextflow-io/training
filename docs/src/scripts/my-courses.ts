// My courses: the learner's picked courses and progress, shared by the header, homepage and
// course pages. Everything is keyed by content slug. Storage sits behind `backend`, so an
// account-backed store can replace localStorage later.

export type Item = { id: string; title: string; parts: string[]; minutes: number; partial?: boolean; sq?: boolean };
export type Goal = 'run' | 'write' | 'custom';
export type State = {
  version: 1;
  goal: Goal;
  /** Picked items, in canonical order. */
  basket: string[];
  /** Done part slugs, per item. Kept when an item leaves the basket. */
  progress: Record<string, string[]>;
};
export type User = { name: string; email: string; initials: string };
export type Backend = { load(): unknown; save(state: State): void };

const KEY = 'nf-training-my-courses';
// TODO(redesign): a stand-in for Seqera login. `?mycourses=loggedin` turns it on, `?mycourses=loggedout` off.
const DEMO_KEY = 'nf-training-my-courses-demo-login';
const DEMO_USER: User = { name: 'Demo learner', email: 'Logged-in preview', initials: 'DL' };
const GOALS: Goal[] = ['run', 'write', 'custom'];

export const localBackend: Backend = {
  load() {
    try {
      return JSON.parse(localStorage.getItem(KEY) ?? 'null');
    } catch {
      return null;
    }
  },
  save(state) {
    try {
      localStorage.setItem(KEY, JSON.stringify(state));
    } catch {}
  },
};

export function formatMinutes(m: number) {
  const h = Math.floor(m / 60);
  const r = m % 60;
  return h && r ? `${h} h ${r} min` : h ? `${h} h` : `${r} min`;
}

/** Drop unknown ids and parts, and put the basket in canonical order. */
export function validate(raw: unknown, items: Item[]): State {
  const s = (raw && typeof raw === 'object' ? raw : {}) as Partial<State>;
  const byId = new Map(items.map((i) => [i.id, i]));
  const basket = Array.isArray(s.basket) ? s.basket : [];
  const progress: State['progress'] = {};
  for (const [id, parts] of Object.entries(s.progress && typeof s.progress === 'object' ? s.progress : {})) {
    const item = byId.get(id);
    const done = Array.isArray(parts) && item ? item.parts.filter((p) => parts.includes(p)) : [];
    if (done.length) progress[id] = done;
  }
  return {
    version: 1,
    goal: GOALS.includes(s.goal as Goal) ? (s.goal as Goal) : 'run',
    basket: items.map((i) => i.id).filter((id) => basket.includes(id)),
    progress,
  };
}

export type Change = { type: 'add' | 'remove' | 'progress' | 'goal' | 'reset' | 'user' | 'sync'; id?: string };

export function createStore(items: Item[], backend: Backend = localBackend) {
  const byId = new Map(items.map((i) => [i.id, i]));
  const order = items.map((i) => i.id);
  let state = validate(backend.load(), items);
  let user = readUser();
  const listeners = new Set<(c: Change) => void>();

  function readUser(): User | null {
    try {
      const flag = new URLSearchParams(location.search).get('mycourses');
      if (flag === 'loggedin') localStorage.setItem(DEMO_KEY, '1');
      if (flag === 'loggedout') localStorage.removeItem(DEMO_KEY);
      return localStorage.getItem(DEMO_KEY) ? DEMO_USER : null;
    } catch {
      return null;
    }
  }
  function commit(change: Change) {
    backend.save(state);
    listeners.forEach((fn) => fn(change));
  }
  const has = (id: string) => state.basket.includes(id);
  const sortBasket = (ids: string[]) => order.filter((id) => ids.includes(id));

  const store = {
    items,
    get: (id: string) => byId.get(id),
    get state() {
      return state;
    },
    get user() {
      return user;
    },
    has,
    basket: () => [...state.basket],
    done: (id: string) => state.progress[id]?.length ?? 0,
    isPartDone: (id: string, part: string) => !!state.progress[id]?.includes(part),
    fraction(id: string) {
      const n = byId.get(id)?.parts.length ?? 0;
      return n ? store.done(id) / n : 0;
    },
    /** Items with progress that are not in the basket. */
    hidden: () => order.filter((id) => !has(id) && store.done(id) > 0),
    overall() {
      let d = 0;
      let p = 0;
      for (const id of state.basket) {
        d += store.done(id);
        p += byId.get(id)!.parts.length;
      }
      return p ? d / p : 0;
    },
    subscribe(fn: (c: Change) => void) {
      listeners.add(fn);
      return () => listeners.delete(fn);
    },
    add(id: string) {
      if (!byId.has(id) || has(id)) return;
      state = { ...state, basket: sortBasket([...state.basket, id]) };
      commit({ type: 'add', id });
    },
    remove(id: string) {
      if (!has(id)) return;
      state = { ...state, basket: state.basket.filter((x) => x !== id) };
      commit({ type: 'remove', id });
    },
    toggle(id: string) {
      has(id) ? store.remove(id) : store.add(id);
    },
    setGoal(goal: Goal) {
      if (goal === state.goal || !GOALS.includes(goal)) return;
      state = { ...state, goal };
      commit({ type: 'goal' });
    },
    /** Mark parts done or not done. Marking anything done adds the item to My courses. */
    setParts(id: string, parts: string[], done: boolean) {
      const item = byId.get(id);
      if (!item) return;
      const now = new Set(state.progress[id] ?? []);
      parts.forEach((p) => (done ? now.add(p) : now.delete(p)));
      const list = item.parts.filter((p) => now.has(p));
      const progress = { ...state.progress };
      if (list.length) progress[id] = list;
      else delete progress[id];
      const added = done && !has(id);
      state = { ...state, progress, basket: added ? sortBasket([...state.basket, id]) : state.basket };
      backend.save(state);
      if (added) listeners.forEach((fn) => fn({ type: 'add', id }));
      listeners.forEach((fn) => fn({ type: 'progress', id }));
    },
    /** Empty the list. Progress is kept. */
    reset() {
      state = { ...state, basket: [] };
      commit({ type: 'reset' });
    },
    login() {
      try {
        localStorage.setItem(DEMO_KEY, '1');
      } catch {}
      user = DEMO_USER;
      listeners.forEach((fn) => fn({ type: 'user' }));
    },
    logout() {
      try {
        localStorage.removeItem(DEMO_KEY);
      } catch {}
      user = null;
      listeners.forEach((fn) => fn({ type: 'user' }));
    },
  };

  // Another tab changed the store.
  globalThis.addEventListener?.('storage', (e: StorageEvent) => {
    if (e.key !== KEY && e.key !== DEMO_KEY) return;
    state = validate(backend.load(), items);
    user = readUser();
    listeners.forEach((fn) => fn({ type: 'sync' }));
  });
  return store;
}

export type Store = ReturnType<typeof createStore>;

let shared: Store | undefined;

/** The page's store, built from the catalogue the header embeds. */
export function myCourses(): Store {
  if (shared) return shared;
  const el = document.getElementById('nf-catalog');
  const items: Item[] = el ? JSON.parse(el.textContent || '[]') : [];
  return (shared = createStore(items));
}
