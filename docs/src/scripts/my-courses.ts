// My courses: the learner's picked courses and progress, shared by the header, homepage and
// course pages. Everything is keyed by content slug. The browser copy lives in localStorage;
// when the learner is logged in with Seqera it mirrors their account through the Router Worker's /api.

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
export type Account = { sub: string; name: string; email: string };
/** What localStorage holds: the state, plus which account it mirrors and whether it has unsaved changes. */
type Stored = State & { account?: Account; pending?: boolean };
export type Backend = { load(): unknown; save(state: Stored): void; clear(): void };
export type Remote = {
  /** `null` when logged out. Throws when offline. */
  load(): Promise<{ user: Account; state: unknown } | null>;
  /** `false` when logged out. Throws when offline. */
  save(state: State): Promise<boolean>;
  logout(): void;
};

const KEY = 'nf-training-my-courses';
const GOALS: Goal[] = ['run', 'write', 'custom'];
const PUSH_DELAY = 1000;

/** Only deploy builds have the Worker's /api; PR previews and plain `astro dev` don't. */
export const loginEnabled = import.meta.env?.PUBLIC_LOGIN === 'true';

export function loginHref() {
  return `/api/auth/login?return=${encodeURIComponent(location.pathname + location.search)}`;
}

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
  clear() {
    try {
      localStorage.removeItem(KEY);
    } catch {}
  },
};

export const accountRemote: Remote = {
  async load() {
    const res = await fetch('/api/progress', { credentials: 'same-origin' });
    if (res.status === 401) return null;
    if (!res.ok) throw new Error(`GET /api/progress: ${res.status}`);
    return res.json();
  },
  async save(state) {
    const res = await fetch('/api/progress', {
      method: 'PUT',
      credentials: 'same-origin',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(state),
      keepalive: true,
    });
    if (res.status === 401) return false;
    if (!res.ok) throw new Error(`PUT /api/progress: ${res.status}`);
    return true;
  },
  logout() {
    location.assign(`/api/auth/logout?return=${encodeURIComponent(location.pathname + location.search)}`);
  },
};

function loggedInHere() {
  try {
    return loginEnabled && document.cookie.split('; ').includes('nf_login=1');
  } catch {
    return false;
  }
}

export function formatMinutes(m: number) {
  const h = Math.floor(m / 60);
  const r = m % 60;
  return h && r ? `${h} h ${r} min` : h ? `${h} h` : `${r} min`;
}

const strings = (x: unknown): string[] => (Array.isArray(x) ? x.filter((s) => typeof s === 'string') : []);

/** Shape only: keeps ids and parts this version's catalogue doesn't know. */
function parse(raw: unknown): State {
  const s = (raw && typeof raw === 'object' ? raw : {}) as Partial<State>;
  const progress: State['progress'] = {};
  for (const [id, parts] of Object.entries(s.progress && typeof s.progress === 'object' ? s.progress : {})) {
    const list = strings(parts);
    if (list.length) progress[id] = list;
  }
  return {
    version: 1,
    goal: GOALS.includes(s.goal as Goal) ? (s.goal as Goal) : 'run',
    basket: strings(s.basket),
    progress,
  };
}

/** Drop unknown ids and parts, and put the basket in canonical order. */
export function validate(raw: unknown, items: Item[]): State {
  const s = parse(raw);
  const byId = new Map(items.map((i) => [i.id, i]));
  const progress: State['progress'] = {};
  for (const [id, parts] of Object.entries(s.progress)) {
    const done = byId.get(id)?.parts.filter((p) => parts.includes(p)) ?? [];
    if (done.length) progress[id] = done;
  }
  return { ...s, basket: items.map((i) => i.id).filter((id) => s.basket.includes(id)), progress };
}

/**
 * `state` plus whatever in `raw` this version's catalogue doesn't know. Every version of the site
 * shares one store, so an older version must not erase progress on courses and parts it lacks.
 */
export function preserve(state: State, raw: unknown, items: Item[]): State {
  const r = parse(raw);
  const known = new Map(items.map((i) => [i.id, new Set(i.parts)]));
  const progress = { ...state.progress };
  for (const [id, parts] of Object.entries(r.progress)) {
    const extra = parts.filter((p) => !known.get(id)?.has(p));
    if (extra.length) progress[id] = [...(progress[id] ?? []), ...extra];
  }
  return { ...state, basket: [...state.basket, ...r.basket.filter((id) => !known.has(id))], progress };
}

/** Union of two states: picked courses and done parts from both. The account's goal wins. */
export function merge(account: unknown, browser: unknown): State {
  const a = parse(account);
  const b = parse(browser);
  const progress = { ...a.progress };
  for (const [id, parts] of Object.entries(b.progress))
    progress[id] = [...new Set([...(progress[id] ?? []), ...parts])];
  return { ...a, goal: account ? a.goal : b.goal, basket: [...new Set([...a.basket, ...b.basket])], progress };
}

const initials = (name: string) =>
  name
    .split(/\s+/)
    .filter(Boolean)
    .map((w) => w[0])
    .filter((_, i, all) => i === 0 || i === all.length - 1)
    .join('')
    .toUpperCase() || '?';
const toUser = (a: Account): User => ({
  name: a.name || a.email,
  email: a.email,
  initials: initials(a.name || a.email),
});

export type Change = { type: 'add' | 'remove' | 'progress' | 'goal' | 'reset' | 'user' | 'sync'; id?: string };

export function createStore(
  items: Item[],
  backend: Backend = localBackend,
  remote: Remote | null = loggedInHere() ? accountRemote : null,
) {
  const byId = new Map(items.map((i) => [i.id, i]));
  const order = items.map((i) => i.id);
  let stored = (backend.load() ?? {}) as Stored;
  // The session has expired: its copy of the account goes, as on logout.
  if (loginEnabled && !remote && stored.account) {
    backend.clear();
    stored = {} as Stored;
  }
  let raw = parse(stored);
  let state = validate(raw, items);
  let account = stored.account;
  let pending = !!stored.pending;
  let user = remote && account ? toUser(account) : null;
  let synced = false;
  let timer: ReturnType<typeof setTimeout> | undefined;
  const listeners = new Set<(c: Change) => void>();
  const emit = (c: Change) => listeners.forEach((fn) => fn(c));

  function save() {
    raw = preserve(state, raw, items);
    backend.save({ ...raw, ...(account && { account }), ...(pending && { pending }) });
  }
  function commit(change: Change) {
    changed();
    emit(change);
  }
  function changed() {
    if (remote || account) {
      pending = true;
      if (synced) {
        clearTimeout(timer);
        timer = setTimeout(push, PUSH_DELAY);
      }
    }
    save();
  }
  async function push() {
    clearTimeout(timer);
    timer = undefined;
    if (!remote || !synced || !pending) return;
    try {
      if (!(await remote.save(raw))) return signedOut();
      pending = false;
      save();
    } catch {}
  }
  /** The session is gone. A copy of the account goes with it, as on logout; anonymous progress stays. */
  function signedOut() {
    remote = null;
    synced = false;
    if (account) {
      backend.clear();
      raw = parse(null);
      state = validate(raw, items);
    }
    account = undefined;
    pending = false;
    user = null;
    emit({ type: 'user' });
  }

  /**
   * Load the account. A browser copy that doesn't mirror this account, or has changes the account
   * hasn't seen, is merged into it once; otherwise the account replaces it, so a part un-marked on
   * another device stays un-marked here.
   */
  async function sync() {
    let res: Awaited<ReturnType<Remote['load']>>;
    try {
      res = await remote!.load();
    } catch {
      return;
    }
    if (!res) return signedOut();
    const foreign = !!account && account.sub !== res.user.sub;
    const mine = !!account && !foreign && !pending && res.state;
    raw = foreign ? parse(res.state) : mine ? parse(res.state) : merge(res.state, raw);
    state = validate(raw, items);
    account = res.user;
    user = toUser(res.user);
    pending = !foreign && !mine;
    synced = true;
    save();
    emit({ type: 'sync' });
    if (pending) await push();
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
    /** Settles once the account has loaded (at once when logged out). */
    ready: remote ? sync() : Promise.resolve(),
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
      changed();
      if (added) emit({ type: 'add', id });
      emit({ type: 'progress', id });
    },
    /** Empty the list. Progress is kept. */
    reset() {
      state = { ...state, basket: [] };
      commit({ type: 'reset' });
    },
    /** Send unsaved changes to the account now. */
    flush: () => push(),
    /** Logging out leaves nothing behind in this browser: the account keeps the learner's courses. */
    async logout() {
      await push();
      const r = remote;
      remote = null;
      backend.clear();
      raw = parse(null);
      state = validate(raw, items);
      account = undefined;
      pending = false;
      user = null;
      emit({ type: 'user' });
      r?.logout();
    },
  };

  // Another tab changed the store.
  globalThis.addEventListener?.('storage', (e: StorageEvent) => {
    if (e.key !== KEY) return;
    stored = (backend.load() ?? {}) as Stored;
    raw = parse(stored);
    state = validate(raw, items);
    pending = !!stored.pending;
    account = stored.account;
    user = remote && account ? toUser(account) : null;
    emit({ type: 'sync' });
  });
  globalThis.addEventListener?.('pagehide', () => timer && push());
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
