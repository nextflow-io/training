# Acceptable vs. Flaggable Output Differences

When comparing tutorial output with documented examples, use these guidelines.

## Acceptable Differences (do not flag)

These vary between runs and are expected:

- **Work directory hashes vs your actual run**: When you run a command, your hash will differ from the documented hash - that's fine
- **Run names**: e.g., `[goofy_torvalds]`, `[happy_turing]`
- **Timestamps**: dates, times, durations
- **Memory/resource values**: may vary by system
- **Nextflow version**: minor version differences (e.g., `25.04.0` vs `25.04.1`)
- **Process ordering**: parallel processes may complete in different order
- **Progress percentages**: timing-dependent

## Flag as Issues

These indicate documentation or code problems:

- **Different text content**: wrong messages, labels, or descriptions
- **Missing or extra output lines**: structural differences
- **Failed commands**: non-zero exit codes when success expected
- **Missing expected files**: outputs not created
- **Wrong file contents**: data doesn't match documentation
- **Different process names**: indicates code/doc mismatch
- **Error messages**: any unexpected errors or warnings
- **Different channel cardinality**: wrong number of items emitted

## Hash Consistency Within Documentation (CRITICAL)

While your actual run hashes will differ from documented hashes (acceptable), the documentation itself must be internally consistent:

### Must Flag These Hash Issues

1. **Invalid hexadecimal characters**: Hashes like `[j6/cdfa66]` or `[5e/4358gc]` are invalid - only `0-9` and `a-f` are valid hex
2. **Duplicate hashes for different non-cached runs**: If the docs show multiple `nextflow run` commands without `-resume`, each run MUST have unique hashes
3. **Prose references not matching terminal output**: If text says "look at hash `[a3/7be2fa]`", that hash must appear in the nearby terminal output
4. **Wrong hashes after `-resume`**: When using `-resume`, cached processes should show the SAME hash as the original run (with "cached: N" status)

### Hash Behavior Rules

| Scenario | Hash Behavior |
|----------|---------------|
| Fresh run (no `-resume`) | New unique hash for every process |
| `-resume` with cache hit | Same hash as original, shows "cached: N" |
| `-resume` with cache miss | New hash (inputs changed) |
| Same run shown twice in docs | Same hash (it's the same execution) |
| Different runs in docs | Different hashes (different executions) |

### How to Verify

```bash
# Find all hashes in a file
grep -oE '\[[a-z0-9]{2}/[a-z0-9]{6}\]' file.md | sort | uniq -c | sort -rn

# Check for invalid hex characters
grep -oE '\[[^]]+\]' file.md | grep -E '\[[^0-9a-f/\]]'
```

Duplicate hash counts > 1 are only acceptable when:
- Same terminal output is shown multiple times
- Prose references the same hash from nearby output
- `-resume` shows cached results from a previous run in the same lesson

## Edge Cases

### Version-dependent output

Some output changes between Nextflow versions. If the documented version differs from the test version:
1. Note the version difference
2. Flag as "version-dependent" rather than "error"
3. Recommend updating documentation if using newer version

### Platform-dependent output

ARM vs x86, macOS vs Linux may show:
- Different container pull messages
- Different resource allocations
- Platform-specific paths

Flag these as "platform-dependent" for review.
