import { defineRouteMiddleware } from '@astrojs/starlight/route-data';
import { pageContext } from './catalog';

/** Each language's homepage is the custom homepage: no table of contents, sidebar on phones only. */
export const isHomepage = (route: { id: string; locale?: string | undefined }) =>
  route.id === '' || route.id === route.locale;

export const onRequest = defineRouteMiddleware(async ({ locals }) => {
  const route = locals.starlightRoute;
  if (isHomepage(route)) {
    // The sidebar stays for the phone menu; Homepage.astro hides it on wider screens.
    route.toc = undefined;
    // A hero makes Starlight skip the page title; the Hero override renders the homepage.
    route.entry.data.hero = { actions: [] };
    return;
  }

  const ctx = await pageContext(route);
  const { course } = ctx;
  if (!course) return;
  const t = locals.t;

  // Pagination stays inside the course: no "previous" before the overview, nothing after the summary,
  // and the feedback survey (linked from the summary) is skipped.
  const pg = route.pagination;
  const isSurvey = (link?: { href: string }) => !!course.survey && !!link?.href.endsWith(`/${course.survey}/`);
  const last = course.parts.at(-1);
  if (ctx.role === 'overview') pg.prev = undefined;
  if (ctx.role === 'summary') pg.next = undefined;
  if (isSurvey(pg.next) && course.summary)
    pg.next = { ...pg.next!, href: ctx.href(course.summary), label: t('nf.sb.summary') };
  if (isSurvey(pg.prev) && last) pg.prev = { ...pg.prev!, href: ctx.href(last.slug), label: last.pageTitle };

  // Overview pages: the table of contents follows the sections the overview renders.
  if (ctx.role !== 'overview') return;
  const toc = (route.toc ??= {
    minHeadingLevel: 2,
    maxHeadingLevel: 3,
    items: [
      { depth: 2, slug: '_top', text: t('tableOfContents.overview'), children: [] },
      ...route.headings
        .filter((h) => h.depth === 2)
        .map((h) => ({ depth: 2, slug: h.slug, text: h.text, children: [] })),
    ],
  });
  const info = (route.entry.data as { additional_information?: Record<string, unknown[]> }).additional_information;
  const item = (slug: string, text: string) => ({ depth: 2, slug, text, children: [] });
  const [top, ...rest] = toc.items;
  const body = rest
    .filter((i) => !['lesson-plan', 'what-youll-learn', 'before-you-start'].includes(i.slug))
    .map((i) => ({ ...i, children: i.children.filter((c) => c.slug !== 'lesson-plan') }));
  toc.items = [
    top,
    ...(info?.learning_objectives?.length ? [item('what-youll-learn', t('nf.course.learn'))] : []),
    ...(info?.audience_prerequisites?.length ? [item('before-you-start', t('nf.course.before'))] : []),
    ...body,
    item('lesson-plan-parts', t('nf.course.plan')),
    ...(course.next.length ? [item('after-this-course', t('nf.course.after'))] : []),
  ];
});
