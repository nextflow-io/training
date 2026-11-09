---
name: add-exercise
description: Add an interactive exercise with solution to an existing lesson. Creates properly formatted exercise and solution collapsible blocks. Use when a lesson needs hands-on practice.
---

Add an exercise with solution to an existing training lesson.

See [../shared/repo-conventions.md](../shared/repo-conventions.md) for directory structure and file conventions.

Follow these steps:

1. Ask the user:

   - Which lesson file to add the exercise to
   - What the exercise should teach/practice
   - Difficulty level (beginner, intermediate, advanced)

2. Read the lesson file to understand context

3. Create the exercise as a collapsible `<details class="exercise">` block, with the solution nested inside it as `<details class="solution">`.
   Leave a blank line after each `<summary>` line and before each `</details>`, or the body does not render as Markdown:

   ````mdx
   <details class="exercise">
   <summary>Exercise: [Title]</summary>

   Brief description of what to do.

   1. First step
   2. Second step
   3. Expected outcome

   :::tip
   Helpful hint if needed
   :::

   <details class="solution">
   <summary>Solution</summary>

   Explanation of the solution approach.

   ```nextflow title="filename.nf" showLineNumbers
   [solution code]
   ```

   ```bash
   nextflow run filename.nf
   # [!output]
   [expected output]
   ```

   </details>

   </details>
   ````

   Add `open` (`<details class="exercise" open>`) only if the exercise should show without a click.

4. If the exercise requires a solution file:

   - Create the solution script in `[module]/solutions/`
   - Follow naming convention matching the exercise
   - Ensure the solution script is complete and tested
   - Add comments explaining key parts

5. Insert the exercise at an appropriate point in the lesson:

   - After explaining the relevant concept
   - Before the section's Takeaway
   - With clear connection to learning objectives

6. Test the exercise:

   - Verify instructions are clear
   - Confirm solution code works
   - Check that exercise difficulty is appropriate
   - Ensure expected output is correct

7. Remind the user to:
   - Test the solution code themselves
   - Consider if hints are needed
   - Preview how it renders in the docs (`/preview`)
   - Check that file paths in exercise are correct
