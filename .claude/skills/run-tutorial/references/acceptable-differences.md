# Acceptable vs. Flaggable Output Differences

When comparing tutorial output with documented examples, use these guidelines.

## Acceptable Differences (do not flag)

These vary between runs and are expected:

- **Work directory hashes**: e.g., `work/a1/b2c3d4...` paths
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
