---
description: Test a Nextflow example script and verify it matches documentation
---

Test a Nextflow example script to ensure it works as documented. Follow these steps:

1. Ask the user which Nextflow script to test (or find recently modified .nf files)

2. For the selected script:
   - Read the script and understand what it does
   - Check if there's corresponding documentation in docs/
   - Look for any parameters defined with `params.`
   - Identify expected outputs

3. Run the script:
   ```bash
   cd [script_directory]
   nextflow run [script.nf]
   ```

4. Verify the output:
   - Check the console output matches documentation examples
   - Verify output files are created as expected
   - Check the work directory structure
   - If publishDir is used, verify results/ directory contents

5. Test with resume:
   ```bash
   nextflow run [script.nf] -resume
   ```
   - Verify it shows "cached" for processes

6. If the script uses parameters, test with different values:
   ```bash
   nextflow run [script.nf] --param value
   ```

7. Compare actual output with documentation:
   - Find references to this script in markdown files
   - Check that documented output matches actual output
   - Check that command examples are correct
   - Verify console output examples are accurate

8. Report:
   - Whether the script runs successfully
   - Any discrepancies between docs and actual behavior
   - Suggestions for documentation improvements
   - Any issues that need fixing
