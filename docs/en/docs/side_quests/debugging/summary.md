# Summary

You have completed the Troubleshooting Workflows mini-course.
This page recaps what you learned in each part and serves as a quick reference for future debugging sessions.

---

## What you learned

### Part 1: Common errors and how to fix them

You worked through the most frequent Nextflow error categories, learning the shape of each error message and the typical fix.
You covered syntax errors (missing braces, wrong keywords, undefined variables, Bash vs Groovy variable handling, statements outside the workflow block), channel structure errors, and process structure errors (missing outputs, missing software, bad resource configuration).

### Part 2: The Nextflow debugging toolkit

Using a single spine pipeline you applied six techniques in sequence:

- **Work-directory forensics** — `.command.sh`, `.command.err`, `.command.out`, `.exitcode` and a plain `ls` to reconstruct exactly what a failed process did
- **`-preview`** — fast parse-time validation that catches syntax errors and bad workflow structure before anything expensive runs
- **`debug true`** — streaming process output and echoes to your terminal as the process runs
- **`-stub-run`** — using stub directives to iterate on workflow logic without the real commands or their containers
- **`-dump-hashes`** — finding out exactly which input component changed when `-resume` re-runs more than you expected
- **A four-phase systematic method** that combines all of the above

---

## Quick reference

### Error shapes (Part 1)

| Hallmark                                              | Likely cause                                |
| ----------------------------------------------------- | ------------------------------------------- |
| `Unexpected input: '<EOF>'`                           | Missing closing brace                       |
| `Invalid process definition`                          | Wrong section keyword (`inputs` vs `input`) |
| `` `name` is not defined ``                           | Typo or missing Groovy variable             |
| `is not defined` on a name set in Bash                | Bash variable not escaped (`\$name`)        |
| `Statements cannot be mixed with script declarations` | Channel definition outside `workflow { }`   |
| `Process requires N channels, M were specified`       | Cardinality mismatch at process call        |
| Stray `[ ]` or `Path value cannot be null`            | Channel emits wrong tuple shape             |
| `Missing output file(s)`                              | Output declaration vs script filename       |
| Exit code 127, `command not found`                    | Software not in container/conda             |
| Exit code 137, `exceeded running time limit`          | Resource directive too tight                |

### Tool cheat sheet (Part 2)

| Tool                     | Reach for it when…                                                                         |
| ------------------------ | ------------------------------------------------------------------------------------------ |
| `work/<hash>/.command.*` | A process failed and you need to see exactly what ran                                      |
| `-preview`               | You just edited the workflow and want a fast structural sanity check                       |
| `debug true` + `echo`    | A process completes but produces something unexpected                                      |
| `-stub-run`              | You're iterating downstream and don't want to wait on slow or containerised upstream tasks |
| `-dump-hashes`           | `-resume` re-ran tasks you thought were cached                                             |
| `nextflow.log`           | Anything that doesn't show up in the terminal output                                       |

### The four-phase method

1. **Parse first** — `nextflow run workflow.nf -preview`
2. **Read the error** — categorise as structural, runtime, or resource
3. **Investigate** — walk the work directory, add `.view()` / `debug true`, use `-stub-run` if iterating
4. **Fix and verify** — minimal change, re-run with `-resume`, use `-dump-hashes` if cache misses surprise you

---

## What's most worth remembering

Three things that are easy to forget:

1. **The work directory has everything.** Any time a process fails, your first move is to find the work directory in the error message and `ls` it.
2. **Any text change inside `script:` busts the cache, including comments and whitespace.** Resource directives don't.
3. **A change upstream can cascade into cache misses downstream**, even when the downstream process is untouched, because its inputs are now different files.
