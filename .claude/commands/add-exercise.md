---
description: Add an exercise with solution to an existing lesson
---

Add an exercise with solution to an existing training lesson. Follow these steps:

1. Ask the user:

   - Which lesson file to add the exercise to
   - What the exercise should teach/practice
   - Difficulty level (beginner, intermediate, advanced)

2. Read the lesson file to understand context

3. Create the exercise using the `??? exercise` admonition:

   ```markdown
   ??? exercise "Exercise: [Title]"

       Brief description of what to do.

       1. First step
       2. Second step
       3. Expected outcome

       !!! tip

           Helpful hint if needed
   ```

4. Create the corresponding solution using `??? solution` or `??? result`:

   ````markdown
   ??? solution "Solution"

       Explanation of the solution approach.

       ```groovy title="filename.nf" linenums="1"
       [solution code]
       ```

       Expected output:

       ```console title="Output"
       [expected output]
       ```
   ````

5. If the exercise requires a solution file:

   - Create the solution script in `[module]/solutions/`
   - Follow naming convention matching the exercise
   - Ensure the solution script is complete and tested
   - Add comments explaining key parts

6. Insert the exercise at an appropriate point in the lesson:

   - After explaining the relevant concept
   - Before the section's Takeaway
   - With clear connection to learning objectives

7. Test the exercise:

   - Verify instructions are clear
   - Confirm solution code works
   - Check that exercise difficulty is appropriate
   - Ensure expected output is correct

8. Remind the user to:
   - Test the solution code themselves
   - Consider if hints are needed
   - Preview how it renders in the docs
   - Check that file paths in exercise are correct
