"""Converter rules, in the order they run. Each module has `apply(page) -> None`.

The order matters:
- frontmatter reads the H1 and intro before anything else rewrites them;
- command_output merges command + output fences before fences rewrites fence lines;
- fences consumes annotation lists before tabs moves code into walkthroughs;
- admonitions turns the remaining indented containers into asides and <details>;
- mdx_escape runs last, over prose only.
"""

from . import admonitions, command_output, components, fences, frontmatter, images, links, mdx_escape, tabs

RULES = [frontmatter, images, command_output, fences, tabs, admonitions, components, links, mdx_escape]
