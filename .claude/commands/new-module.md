---
description: Scaffold a new complete training module with all standard files
---

Create a new training module in this repository. Follow these steps:

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

4. Add the module to mkdocs.yml navigation structure

5. Remind the user to:
   - Create any associated Nextflow example scripts in the root-level directory
   - Add Excalidraw diagrams if needed
   - Test all examples thoroughly
   - Run validation tools before committing
   - Preview the site locally
