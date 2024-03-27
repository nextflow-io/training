# Linting

Linting is a static analysis process that helps ensure code quality by automatically identifying syntax errors, potential bugs, and adherence to coding standards. By enforcing consistency and best practices, linting enhances code readability, reduces errors, and streamlines the development workflow.

It is recommended that you use the `nf-core lint` command to check for inconsistencies in your code, compare your code against source code, and compare your code against nf-core standards.

Executing the `nf-core lint` command from within your pipeline repository will print a list of ignored tests, warnings, failed tests, and a summary.

```console
╭───────────────────────╮
│ LINT RESULTS SUMMARY  │
├───────────────────────┤
│ [✔] 184 Tests Passed  │
│ [?]   0 Tests Ignored │
│ [!]  29 Test Warnings │
│ [✗]   0 Tests Failed  │
╰───────────────────────╯
```

The linting tests in the nf-core template are designed for nf-core pipeline that are shared as a part of the nf-core community. As such, you may find that you want to ignore certain linting failures that are not required for your use case.

A full list of checks, descriptions of how they work, and how they can be customized can be found on the [tools documentation website](https://nf-co.re/tools/docs).

To specify which tests you would like to skip while linting you can modify the `.nf-core.yml` file.

For example, there is a linting test that checks the `CODE_OF_CONDUCT.md` and will throw warnings and errors if it has been edited or removed.

```console
╭─ [✗] 1 Pipeline Test Failed ────────────────────────────────────╮
│                                                                 │
│ files_unchanged: CODE_OF_CONDUCT.md does not match the template │
│                                                                 │
╰─────────────────────────────────────────────────────────────────╯
```

In this scenario, the nf-core lint command will offer a command to fix this warning.

```console
Tip: Some of these linting errors can automatically be resolved with the following command:

    nf-core lint --fix files_unchanged
```

Alternatively, you can add this file to `files_unchanged:` in the `.nf-core.yml` file to ignore the tests:

```yml title=".nf-core.yml"
repository_type: pipeline
lint:
    files_unchanged:
        - CODE_OF_CONDUCT.md
```

If you run nf-core lint again, you will see that the test is now ignored and there are no more failed tests.

```console
╭─ [?] 1 Pipeline Test Ignored ────────────────────────────────────────╮
│                                                                      │
│ files_unchanged: File ignored due to lint config: CODE_OF_CONDUCT.md │
│                                                                      │
╰──────────────────────────────────────────────────────────────────────╯
```

!!! question "Exercise"

    Make an edit to the `CODE_OF_CONDUCT.md` in your pipeline repository (e.g., add another bullet point). Use the `nf-core lint` command to see if it passes or fails. Add the `CODE_OF_CONDUCT.md` to the files unchanged in your `.nf-core.yml` file and lint your pipeline again to show that the test has been ignored.

!!! question "Bonus Exercise"

    Modify the `.nf-core.yml` file to prevent `pipeline_todos` from showing as warnings in your lint tests.

    ??? "Solution"

        ```yml title=".nf-core.yml"
        repository_type: pipeline
        lint:
            pipeline_todos: false
            files_unchanged:
                - CODE_OF_CONDUCT.md
        ```
