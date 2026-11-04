---
name: stop-preview
description: Stop a running Astro dev or preview server. Use when port 4321 is busy, switching tasks, or done editing documentation.
---

Stop the Astro dev server (and any `astro preview` server) that is running in the background.

See [../shared/repo-conventions.md](../shared/repo-conventions.md) for directory structure.

This will:

1. Find any running background dev or preview server
2. Stop it
3. Confirm the server has been stopped

If no preview server is running, you'll be notified.

Use the following steps, in `docs/`:

1. Run `npx astro dev status` and `npx astro preview status` to find running servers
2. Stop each running one with `npx astro dev stop` or `npx astro preview stop`
3. Run the status commands again to confirm they have been stopped

If something still holds port 4321 (a server started in another terminal), find it with `lsof -i :4321` and ask the user before stopping it.
