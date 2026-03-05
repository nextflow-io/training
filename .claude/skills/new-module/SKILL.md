---
name: new-module
description: Scaffold a complete new training module with standard directory structure, index, orientation, lessons, survey, and next-steps files. Use when creating a new course like 'Advanced Pipelines' or a new domain module.
---

Create a new training module in this repository.

See [../shared/repo-conventions.md](../shared/repo-conventions.md) for directory structure and file conventions.

Follow these steps:

1. Ask the user for:

   - Module name and location (under docs/)
   - Module title for display
   - Brief description of what the module covers

2. Create the complete module structure:

   ```
   docs/[module_name]/
   ├── index.md (overview with learning objectives)
   ├── 00_orientation.md (prerequisites, setup, introduction)
   ├── 01_first_topic.md (first lesson)
   ├── survey.md (feedback collection)
   ├── next_steps.md (what to learn after completing module)
   ├── img/ (directory for images and diagrams)
   └── solutions/ (directory for solution code)
   ```

3. Each file should follow repository conventions:

   - Proper heading numbering with trailing periods
   - Takeaway and What's next? sections
   - Appropriate use of admonitions
   - Code examples with proper formatting

4. Add the module to `docs/en/mkdocs.yml` navigation structure

5. Remind the user to:
   - Create any associated Nextflow example scripts in the root-level directory
   - Add Excalidraw diagrams if needed
   - Test all examples thoroughly
   - Run validation tools before committing
   - Preview the site locally
   - Translations are handled automatically - when merged to master, the translation workflow will create PRs for each language
