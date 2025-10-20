---
description: Update mkdocs.yml navigation after adding or moving content
---

Update the navigation structure in mkdocs.yml after adding, removing, or reorganizing training content.

1. First, identify what changed:
   - New lesson/module files added?
   - Files moved or renamed?
   - Files deleted?

2. Read the current mkdocs.yml nav section

3. Show the user where the changes need to be made in the nav structure

4. Explain the nav structure format:
   ```yaml
   nav:
     - Module Name:
       - module/index.md
       - module/00_orientation.md
       - module/01_lesson.md
       - module/02_lesson.md
   ```

5. For new content:
   - Determine the correct location in the hierarchy
   - Maintain numbered sequence
   - Use proper indentation (2 spaces)
   - Use meaningful display names (not just filenames)

6. Update the nav section preserving:
   - Existing structure and organization
   - Proper YAML formatting
   - Alphabetical or logical ordering
   - Consistent indentation

7. Verify:
   - All new files are included
   - No broken references
   - Proper nesting for submodules
   - Display names are clear and helpful

8. After updating, remind user to:
   - Preview the site to verify navigation works
   - Check that new pages are accessible
   - Verify breadcrumbs and next/previous links work correctly
   - Test search functionality includes new content
