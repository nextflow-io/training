# Linting

Linting is a static analysis process that helps ensure code quality by automatically identifying syntax errors, potential bugs, and adherence to coding standards. By enforcing consistency and best practices, linting enhances code readability, reduces errors, and streamlines the development workflow.

As a part of nf-core tools, the `nf-core lint` command can be used to check for inconsistencies in your code, compare your code against source code, and compare your code against nf-core standards.

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

The linting tests in the nf-core template are primarily designed for pipelines that are shared as a part of the nf-core project. As such, you may find that you want to ignore certain linting failures if they are not required.

A full list of checks, descriptions of how they work, and how they can be customized can be found on the [tools documentation website](https://nf-co.re/tools/docs).

The `.nf-core.yml` file is contained within the nf-core pipeline template and can be modified to choose which linting tests you would like to skip.

For example, by default, there is a linting test that check if the `CODE_OF_CONDUCT.md` matches the template. If the `CODE_OF_CONDUCT.md` has been modified the linting test will fail.

```console
╭─ [✗] 1 Pipeline Test Failed ────────────────────────────────────╮
│                                                                 │
│ files_unchanged: CODE_OF_CONDUCT.md does not match the template │
│                                                                 │
╰─────────────────────────────────────────────────────────────────╯
```

In this scenario, the `nf-core lint` command will suggest a command to fix this warning.

```console
Tip: Some of these linting errors can automatically be resolved with the following command:

    nf-core lint --fix files_unchanged
```

However, if you wish to remove or modify this file you would need to ignore this test. To do this, you list the `CODE_OF_CONDUCT.md` as a `files_unchanged:` in the `.nf-core.yml` file.

```yml title=".nf-core.yml"
repository_type: pipeline
lint:
    files_unchanged:
        - CODE_OF_CONDUCT.md
```

If you run `nf-core lint` again, you will see that the test is now ignored and there are no more failed tests.

```console
╭─ [?] 1 Pipeline Test Ignored ────────────────────────────────────────╮
│                                                                      │
│ files_unchanged: File ignored due to lint config: CODE_OF_CONDUCT.md │
│                                                                      │
╰──────────────────────────────────────────────────────────────────────╯
```

!!! question "Exercise"

    Make an edit to the `CODE_OF_CONDUCT.md` in your pipeline repository (e.g., add another bullet point). Use the `nf-core lint` command to see if it passes or fails. Add the `CODE_OF_CONDUCT.md` as a `files_unchanged:` in the `.nf-core.yml` file and lint your pipeline again to show that the test has been ignored.

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
