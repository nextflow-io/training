---
name: new-module
description: Scaffold a complete new training module with standard directory structure, index, orientation, lessons, survey, and next-steps files. Use when creating a new course like 'Advanced Pipelines' or a new domain module.
---

Create a new training module in this repository.

See [../shared/repo-conventions.md](../shared/repo-conventions.md) for directory structure and file conventions.

For a new `nf4_science` domain course, copy the template instead: follow `docs/src/content/docs/nf4_science/_template/_README.md`.

Follow these steps:

1. Ask the user for:

   - Module name and location (under `docs/src/content/docs/`)
   - Module title for display
   - Brief description of what the module covers

2. Create the complete module structure:

   ```
   docs/src/content/docs/[module_name]/
   ├── index.mdx (overview with learning objectives)
   ├── 00_orientation.mdx (prerequisites, setup, introduction)
   ├── 01_first_topic.mdx (first lesson)
   ├── survey.mdx (feedback collection)
   ├── next_steps.mdx (what to learn after completing module)
   └── img/ (directory for images and diagrams)

   [module-name]/ (scripts, at the repository root)
   └── solutions/ (directory for solution code)
   ```

3. Each file should follow repository conventions:

   - `title:` and `description:` frontmatter, no H1
   - Proper heading numbering with trailing periods
   - Takeaway and What's next? sections
   - Appropriate use of asides and `<details>` collapsibles
   - Code examples with proper formatting

4. Write `index.mdx` with the index page template (full spec in CONTRIBUTING.md):

   ```mdx
   ---
   title: Module Title
   description: One-sentence summary of the module.
   tableOfContents: false
   page_type: index_page
   index_type: course
   additional_information:
     technical_requirements: true
     learning_objectives:
       - First objective
     audience_prerequisites:
       - '**Audience:** Who this is for'
   ---

   **One-line summary in bold.**

   Summary paragraph.

   <CourseInfo frontmatter={frontmatter} />

   ## Course overview

   ...

   <LinkButton href="00_orientation.mdx" variant="primary" icon="right-arrow">Get started</LinkButton>
   ```

   `index_type` is `course` or `module`. The build validates these fields against `docs/src/content.config.ts`.

5. Add the module to `docs/src/sidebar.mjs` as a group, following the existing ones:

   ```js
   {
     "label": "Module Title",
     "collapsed": true,
     "items": [
       { "label": "Overview", "slug": "module_name" },
       "module_name/00_orientation",
       "module_name/01_first_topic",
       "module_name/next_steps",
       "module_name/survey"
     ]
   }
   ```

   Existing groups carry a `translations` map for their labels; copy the "Overview" one from a sibling group.

6. Remind the user to:
   - Create any associated Nextflow example scripts in the root-level directory
   - Add Excalidraw diagrams (`img/*.excalidraw.svg`) if needed
   - Test all examples thoroughly
   - Run validation tools before committing (`/validate`)
   - Preview the site locally (`/preview`)
   - Translations are handled by the translation workflow, which maintainers run before a release
