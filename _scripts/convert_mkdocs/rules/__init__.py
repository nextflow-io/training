"""Converter rules, in the order they run. Each module has `apply(page) -> None`.

The order matters:
- frontmatter reads the H1 and intro before anything else rewrites them;
- command_output merges command + output fences before fences rewrites fence lines;
- fences consumes annotation lists before tabs moves code into walkthroughs;
- admonitions turns the remaining indented containers into `:::` directives;
- details turns raw HTML <details> blocks into collapsible directives;
- mdx_escape runs last, over prose only.
"""

from . import admonitions, command_output, components, details, fences, frontmatter, images, links, mdx_escape, tabs

RULES = [frontmatter, images, command_output, fences, tabs, admonitions, details, components, links, mdx_escape]
