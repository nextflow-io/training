// Build-time course catalogue: courses.mjs joined with sidebar.mjs and each page's frontmatter.
import { getCollection } from 'astro:content';
import { courses as courseConfig, groups, sideQuests as sqConfig } from './courses.mjs';
import { sidebar } from './sidebar.mjs';

export type Part = { slug: string; title: string; minutes?: number };
export type Course = {
  id: string;
  title: string;
  description: string;
  group: keyof typeof groups;
  icon: string;
  blurb: string;
  level?: string;
  next: string[];
  start?: string;
  summary?: string;
  parts: Part[];
  minutes: number;
  /** Some parts have no `duration`, so `minutes` is a lower bound. */
  partial: boolean;
};
export type Quest = Part & { parts: Part[] };
export type QuestGroup = { label: string; quests: Quest[] };
export type SideQuests = {
  id: string;
  title: string;
  description: string;
  group: keyof typeof groups;
  icon: string;
  start?: string;
  groups: QuestGroup[];
};
/** What the client store needs about each item a learner can add to My courses. */
export type ClientItem = {
  id: string;
  title: string;
  parts: string[];
  minutes: number;
  partial?: boolean;
  sq?: boolean;
};

type SidebarItem = string | { slug?: string; label?: string; link?: string; items?: SidebarItem[] };

const leaf = (s: string) => s.slice(s.lastIndexOf('/') + 1);
const isPart = (slug: string) => /^\d\d_/.test(leaf(slug)) && !/^00_/.test(leaf(slug));
const slugOf = (item: SidebarItem) => (typeof item === 'string' ? item : item.slug);
const flatten = (items: SidebarItem[]): string[] =>
  items.flatMap((i) => (typeof i === 'string' ? [i] : i.items ? flatten(i.items) : i.slug ? [i.slug] : []));

/** The sidebar group whose items include the overview page `id`. */
function findGroup(items: SidebarItem[], id: string): SidebarItem[] | undefined {
  for (const item of items) {
    if (typeof item === 'string' || !item.items) continue;
    if (item.items.some((i) => slugOf(i) === id)) return item.items;
    const found = findGroup(item.items, id);
    if (found) return found;
  }
}

export const shortTitle = (title: string) => title.replace(/^Part \d+:\s*/i, '');

let cache: Promise<{ courses: Course[]; sideQuests: SideQuests }> | undefined;

export function getCatalog() {
  return (cache ??= build());
}

async function build() {
  const entries = new Map((await getCollection('docs')).map((e) => [e.id || 'index', e]));
  const page = (slug: string) => {
    const entry = entries.get(slug);
    if (!entry) throw new Error(`courses.mjs: no page for "${slug}"`);
    return entry.data as { title: string; description?: string; duration?: number };
  };
  const part = (slug: string): Part => ({ slug, title: shortTitle(page(slug).title), minutes: page(slug).duration });

  const courses: Course[] = courseConfig.map((c) => {
    const items = findGroup(sidebar, c.id);
    if (!items) throw new Error(`courses.mjs: "${c.id}" is not an overview page in sidebar.mjs`);
    const slugs = flatten(items);
    const parts = slugs.filter(isPart).map(part);
    return {
      ...c,
      group: c.group as Course['group'],
      next: c.next ?? [],
      title: page(c.id).title,
      description: page(c.id).description ?? '',
      start: slugs.find((s) => /^0+_|^orientation$/.test(leaf(s))),
      summary: slugs.find((s) => leaf(s) === 'next_steps'),
      parts,
      minutes: parts.reduce((n, p) => n + (p.minutes ?? 0), 0),
      partial: parts.some((p) => !p.minutes),
    };
  });

  const sqItems = findGroup(sidebar, sqConfig.id);
  if (!sqItems) throw new Error(`courses.mjs: "${sqConfig.id}" is not an overview page in sidebar.mjs`);
  const quest = (item: SidebarItem): Quest => {
    const slugs = typeof item === 'string' ? [item] : flatten(item.items ?? []);
    const main = part(slugs[0]);
    const parts = slugs.filter(isPart).map(part);
    return { ...main, parts: parts.length ? parts : [main] };
  };
  const questGroups: QuestGroup[] = sqItems
    .filter((i): i is Exclude<SidebarItem, string> => typeof i !== 'string' && !!i.items)
    .map((g) => ({ label: g.label ?? '', quests: g.items!.map(quest) }));

  const sideQuests: SideQuests = {
    ...sqConfig,
    group: sqConfig.group as SideQuests['group'],
    title: page(sqConfig.id).title,
    description: page(sqConfig.id).description ?? '',
    start: flatten(sqItems).find((s) => leaf(s) === 'orientation'),
    groups: questGroups,
  };
  return { courses, sideQuests };
}

/** Courses then Side Quests, in canonical order: the shape the client store uses. */
export async function clientItems(): Promise<ClientItem[]> {
  const { courses, sideQuests } = await getCatalog();
  return [
    ...courses.map((c) => ({
      id: c.id,
      title: c.title,
      parts: c.parts.map((p) => p.slug),
      minutes: c.minutes,
      ...(c.partial && { partial: true }),
    })),
    ...sideQuests.groups.flatMap((g) =>
      g.quests.map((q) => ({
        id: q.slug,
        title: q.title,
        parts: q.parts.map((p) => p.slug),
        minutes: q.minutes ?? 0,
        sq: true,
      })),
    ),
  ];
}

export { formatMinutes } from './scripts/my-courses';
