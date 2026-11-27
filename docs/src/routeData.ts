import { defineRouteMiddleware } from '@astrojs/starlight/route-data';

/** Each language's homepage is the custom homepage: no sidebar, no table of contents. */
export const isHomepage = (route: { id: string; locale?: string | undefined }) =>
  route.id === '' || route.id === route.locale;

export const onRequest = defineRouteMiddleware(({ locals }) => {
  const route = locals.starlightRoute;
  if (!isHomepage(route)) return;
  route.hasSidebar = false;
  route.toc = undefined;
  // A hero makes Starlight skip the page title; the Hero override renders the homepage.
  route.entry.data.hero = { actions: [] };
});
