# Deep Lesson Review Checklist

Use this checklist when reviewing a specific lesson file.

## Structure

- [ ] `title:` and `description:` frontmatter; no H1 in the body
- [ ] Proper heading numbering with trailing periods
- [ ] Heading levels match numbering depth (## for 1., ### for 1.1.)
- [ ] Each major section has "### Takeaway"
- [ ] Each major section has "### What's next?"
- [ ] Logical flow from simple to complex
- [ ] Clear learning objectives stated or implied

## Formatting

- [ ] Code blocks have proper titles, line numbers (`showLineNumbers`, `startLineNumber`), and highlights (`{...}`)
- [ ] Nextflow code uses ` ```nextflow `; real Groovy stays ` ```groovy `
- [ ] Command output sits in the `bash` block after `# [!output]`, with key lines marked `# [!code highlight]`
- [ ] File paths use proper markdown formatting
- [ ] Before/After comparisons use `<CodeWalkthrough>` with `step="Before"` and `step="After"`
- [ ] Asides (`:::note`, `:::tip`, `:::caution`) and `<details>` collapsibles properly formatted, with blank lines around the body
- [ ] Each sentence on new line (for clean diffs)

## Content Accuracy

- [ ] Technical accuracy of Nextflow code
- [ ] Commands are correct and runnable
- [ ] Parameter syntax correct (-- for pipeline, - for Nextflow)
- [ ] Examples progress logically
- [ ] Common pitfalls addressed
- [ ] Edge cases explained

## Teaching Effectiveness

- [ ] Clear explanations for beginners
- [ ] Concepts introduced before use
- [ ] Examples are relevant and motivating
- [ ] Exercises appropriate for skill level
- [ ] Solutions available for exercises
- [ ] Good use of tips and warnings

## Cross-References

- [ ] Links to related lessons are relative `.mdx` links and work (the build's links validator checks this)
- [ ] References to files/scripts are correct
- [ ] External links are valid
- [ ] Prerequisites clearly stated

## Examples & Code

- [ ] All Nextflow examples are syntactically correct
- [ ] Variable names are clear and consistent
- [ ] Comments explain non-obvious code
- [ ] Examples can be run as shown
- [ ] Output examples match what code produces

## Writing Style

See CLAUDE.md for full guidelines.

- [ ] No LLM-style interjections ("Remember when...?", "Don't worry") — "Let's" is acceptable in tutorial prose
- [ ] No exclamation-based emphasis ("Much faster!", "So powerful!")
- [ ] Em-dash elaborations replaced with periods or semicolons
- [ ] List explanations use colons, not hyphens
- [ ] Single clear explanation per concept (not multiple phrasings)
- [ ] Professional, direct tone throughout
