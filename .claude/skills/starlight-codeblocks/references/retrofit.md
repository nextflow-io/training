# Add features to an existing site

Use this file when the plugin goes into a site that already has pages. Also use it when the author asks you to improve the code blocks of existing pages. The aim is a page that teaches better. The aim is not the most features.

## Before you start

1. Build the site before you install the plugin. Keep the list of warnings, so you can tell new warnings from old ones.
2. Install the plugin, then build again. Some features start on their own (see "Features that start on their own" in SKILL.md). Open the pages with Python, Nextflow, `diff` or terminal blocks, and check that these features help there. Turn off a feature for a block only if it makes that block worse.
3. Agree the scope with the author: one page, one sidebar group or the whole site.

## Work one page at a time

Read the whole page before you change a block. Find out what the page teaches and who reads it. A tutorial, a reference page and a troubleshooting page need different features.

Then go through the code blocks in order, one at a time. For each block, answer these questions:

- What does the reader do with this block? They can read it to understand an idea, copy it into a file, run it in a terminal, or compare it with another block.
- What does the prose around the block say about it? Prose that names lines, values or steps points to a feature.
- Which part of the block matters for this page? Often it is a few lines of a longer file.
- Does the reader need to change a value before the code works?

Choose a feature only when an answer points to one. Most blocks need no change, or one feature. Use two or more features in a block only when each one has its own reason.

## Signals in existing pages

| What you find | Consider |
|---|---|
| Prose such as "line 3 sets the port" or "the `cache` section" | Focus on those lines, or code mentions if the prose names more than one part |
| A numbered list under the block that explains its lines | Footnotes or annotations, with each note on its line |
| Code comments that explain the code to the reader, not to a future maintainer | Inline callouts, annotations or footnotes |
| Two blocks, before and after an edit | One `diff` block, or `[!code ++]` and `[!code --]` |
| A long file where the text is about one part | Focus, hidden lines or an expandable block |
| Imports and set-up at the top that the text never mentions | Hidden lines |
| `YOUR_API_KEY`, `<your-project>` or a note such as "replace X with your value" | Fill-in placeholders |
| `<Tabs>` with one code block in each tab, such as npm, pnpm and Yarn | Code tabs |
| A step-by-step tutorial that shows the same file again and again as it grows | Code walkthrough, or scrollycoding if prose explains each step |
| A "wrong" and a "right" version of the same code | Line states: `[!code error]` and `[!code success]` |
| YAML, Python or Makefile code where the text is about indentation | Visible whitespace |
| Names of functions or options in prose, as plain inline code | Inline code highlighting |
| A link in the prose to the reference of a name in the block | A code link on that name |
| Python, JavaScript or TypeScript that the reader wants to try | Run code or open in playground |

## Rules

- Keep the code the same. A change can add attributes and directives, but the copied text stays the same, unless the author asks for a change to the code.
- Keep the prose true. If a feature does the job of some prose, such as a list of line notes, move the text into the feature. Keep the prose if it adds something.
- Keep a code comment that the reader must have in the copied file. Move a comment into a feature only if it is there to explain the code on the page.
- Use one note style in a block, and use the same feature for the same job on all pages of the site.
- Do not add a feature to show that the plugin can do it. If you are not sure that a feature helps the reader, leave the block as it is.

## After each page

1. Build the site and read every warning from `starlight-codeblocks`.
2. Open the page and check each block that you changed. Check that the copy button copies the same code as before.
3. Tell the author what you changed on the page and why, in one line for each block. Name the blocks that you left as they are, if a feature was a near choice.
