# Testing with nf-test

Being able to systematically test that every part of your workflow is doing what it's supposed to do is critical for reproducibility and long-term maintenance, and can be a huge help during the development process.

Let's take a minute to talk about why testing is so important. If you're developing a workflow, one of the first things you will do is grab some test data that you know is valid and should produce a result. You add the first process to the pipeline and wire it up to your inputs to make it work. Then, to check it's all working, so you run it on the test data. Assuming that works, you move on to the next process and run the test data again. You repeat this process until you have a pipeline that you're happy with.

Then, maybe you add a simple true or false parameter such as `--skip_process`. Now you must run the pipeline twice, once with each parameter to make sure it works as expected. But wait, how do we check if the `--skip_process` actually skips the process? We have to dig into the outputs or check the log files! This is a pain and prone to error.

As you develop your pipeline, it will quickly become so complex that manually testing every iteration is slow and error prone. Furthermore, if you do find an error it will be very difficult to pin down exactly where in your pipeline the error is coming from. This is where testing comes in.

Testing allows you to systematically check that every part of your pipeline is working as expected. The benefits to a developer of well written tests are huge:

- **Confidence**: Because the tests cover the entire pipeline, you can be be confident changing something doesn't affect anything else
- **Trust**: When multiple developers work on the pipeline, they know the other developers haven't broken the pipeline and every component.
- **Transparency**: The tests show where a pipeline is failing and make it easier to track down the problem. They also function as a form of documentation, showing how to run a process or workflow.
- **Speed**: Because the tests are automated, they can be run very quickly and repeatedly. You can iterate quickly with less fear of introducing new bugs.

There are lots of different types of tests we can write:

1. **Module-level tests**: For individual processes
2. **Workflow-level tests**: For a single workflow
3. **Pipeline-level tests**: For the pipeline as a whole
4. **Performance tests**: For the speed and efficiency of the pipeline
5. **Stress tests**: Assessing the pipeline's performance under extreme conditions to determine its limits

Testing individual processes is analogous to unit tests in other languages. Testing the workflow or the entire pipeline is analogous to what's called integration tests in other languages, where we test the interactions of the components.

[**nf-test**](https://www.nf-test.com/) is a tool that allows you to write module, workflow and pipeline level test. In short, it allows you to systematically check every individual part of the pipeline is working as expected, _in isolation_.

### Learning goals

In this side quest, you'll learn to use nf-test to write a workflow-level test for the pipeline as well as module-level tests for the three processes it calls on.

By the end of this side quest, you'll be able to use the following techniques effectively:

- Initialize nf-test in your project
- Generate module-level and workflow-level tests
- Add common types of assertions
- Understand when to use snapshots vs. content assertions
- Run tests for an entire project

These skills will help you implement a comprehensive testing strategy in your pipeline projects, ensuring they are more robust and maintainable.

### Prerequisites

Before taking on this side quest, you should:

- Have completed the [Hello Nextflow](../hello_nextflow/README.md) tutorial or equivalent beginner's course.
- Be comfortable using basic Nextflow concepts and mechanisms (processes, channels, operators, working with files, meta data)

---

## 0. Get started

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Move into the project directory

Let's move into the directory where the files for this tutorial are located.

```bash
cd side-quests/nf-test
```

You can set VSCode to focus on this directory:

```bash
code .
```

#### Review the materials

You'll find a main workflow file and a CSV file called `greetings.csv` that contains the input to the pipeline.

```console title="Directory contents"
.
├── greetings.csv
└── main.nf
```

For a detailed description of the files, see the [warmup from Hello Nextflow](../hello_nextflow/00_orientation.md).

The workflow we'll be testing is a subset of the Hello workflow built in [Hello Workflow](../hello_nextflow/03_hello_workflow.md).

??? example "What does the Hello Nextflow workflow do?"

    If you haven't done the [Hello Nextflow](../hello_nextflow/index.md) training, here's a quick overview of what this simple workflow does.

    The workflow takes a CSV file containing greetings, runs four consecutive transformation steps on them, and outputs a single text file containing an ASCII picture of a fun character saying the greetings.

    The four steps are implemented as Nextflow processes (`sayHello`, `convertToUpper`, `collectGreetings`, and `cowpy`) stored in separate module files.

    1. **`sayHello`:** Writes each greeting to its own output file (e.g., "Hello-output.txt")
    2. **`convertToUpper`:** Converts each greeting to uppercase (e.g., "HELLO")
    3. **`collectGreetings`:** Collects all uppercase greetings into a single batch file
    4. **`cowpy`:** Generates ASCII art using the `cowpy` tool

    The results are published to a directory called `results/`, and the final output of the pipeline (when run with default parameters) is a plain text file containing ASCII art of a character saying the uppercased greetings.

    In this side quest, we use an intermediate form of the Hello workflow that only contains the first two processes. <!-- TODO: change this to use the full finished workflow as suggested in https://github.com/nextflow-io/training/issues/735 -->

The subset we'll be working with is composed of two processes: `sayHello` and `convertToUpper`.
You can see the full workflow code below.

??? example "Workflow code"

    ```groovy title="main.nf"
    /*
    * Pipeline parameters
    */
    params.input_file = "greetings.csv"

    /*
    * Use echo to print 'Hello World!' to standard out
    */
    process sayHello {

        publishDir 'results', mode: 'copy'

        input:
            val greeting

        output:
            path "${greeting}-output.txt"

        script:
        """
        echo '$greeting' > '$greeting-output.txt'
        """
    }

    /*
    * Use a text replace utility to convert the greeting to uppercase
    */
    process convertToUpper {

        publishDir 'results', mode: 'copy'

        input:
            path input_file

        output:
            path "UPPER-${input_file}"

        script:
        """
        cat '$input_file' | tr '[a-z]' '[A-Z]' > UPPER-${input_file}
        """
    }

    workflow {

        // create a channel for inputs from a CSV file
        greeting_ch = channel.fromPath(params.input_file).splitCsv().flatten()

        // emit a greeting
        sayHello(greeting_ch)

        // convert the greeting to uppercase
        convertToUpper(sayHello.out)
    }
    ```

#### Run the workflow

Let's run the workflow to make sure it's working as expected.

```bash
nextflow run main.nf
```

```console title="Result of running the workflow"
 N E X T F L O W   ~  version 24.10.2

Launching `main.nf` [soggy_linnaeus] DSL2 - revision: bbf79d5c31

executor >  local (6)
[f7/c3be66] sayHello (3)       | 3 of 3 ✔
[cd/e15303] convertToUpper (3) | 3 of 3 ✔
```

CONGRATULATIONS! You just ran a test!

"Wait, what? I just ran the workflow and it worked! How is that a test?"

Good question!

Let's break down what just happened.

You ran the workflow with the default parameters, you confirmed it worked and you're happy with the results. This is the essence of testing. If you worked through the Hello Nextflow training course, you'll have noticed we always started every section by running the workflow we were using as a starting point, to confirm everything is set up correctly.

Testing software essentially does this process for us.

#### Review the assignment

Your challenge is to add standardized tests to this workflow using nf-test, in order to make it easy to verify that every part continues to work as expected in case any further changes are made.

<!-- TODO: give a bit more details, similar to how it's done in the Metadata side quest -->

#### Readiness checklist

Think you're ready to dive in?

- [ ] I understand the goal of this course and its prerequisites
- [ ] My codespace is up and running
- [ ] I've set my working directory appropriately
- [ ] I've run the workflow successfully
- [ ] I understand the assignment

If you can check all the boxes, you're good to go.

---

## 1. Initialize `nf-test`

The `nf-test` package provides an initialization command that sets up a few things in order for us to start developing tests for our project.

```bash
nf-test init
```

This should produce the following output:

```bash
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

It also creates a `tests` directory containing a configuration file stub.

### 1.1. Generate an nf-test

`nf-test` comes with a set of tools for building nf-test files, saving us the majority of the work. These come under the subcommand `generate`. Let's generate a test for the pipeline:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote pipeline test file '/workspaces/training/side-quests/nf-test/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

This will create a `main.nf.test` file within the `tests` directory. This is our pipeline level test file. If you run `tree tests/` you should see something like this:

```console title="Test directory contents"
tests/
├── main.nf.test
└── nextflow.config
```

The `main.nf.test` file is our pipeline level test file. Let's open it up and take a look at the contents.

```groovy title="tests/main.nf.test"
nextflow_pipeline {

    name "Test Workflow main.nf"
    script "main.nf"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
        }

        then {
            assert workflow.success
        }

    }

}
```

We'll take a second to understand the structure of the test file.

The `nextflow_pipeline` block is the entry point for all pipeline level tests. It contains the following:

- `name`: The name of the test.
- `script`: The path to the pipeline script.

The `test` block is the actual test. It contains the following:

- `when`: The conditions under which the test should be run. This includes the parameters that will be used to run the pipeline.
- `then`: The assertions that should be made. This includes the expected outcomes of the pipeline.

In plain English, the logic of the test reads as follows:
"**When** these _parameters_ are provided to this _pipeline_, **then** we expect to see these results."

This isn't a functional test, we will demonstrate how to turn it into one in the next section.

### A Note on Test Names

In the example above, we used the default name "Should run without failures" which is appropriate for a basic test that just checks if the pipeline runs successfully. However, as we add more specific test cases, we should use more descriptive names that indicate what we're actually testing. For example:

- "Should convert input to uppercase" - when testing specific functionality
- "Should handle empty input gracefully" - when testing edge cases
- "Should respect max memory parameter" - when testing resource constraints
- "Should create expected output files" - when testing file generation

Good test names should:

1. Start with "Should" to make it clear what the expected behavior is
2. Describe the specific functionality or scenario being tested
3. Be clear enough that if the test fails, you know what functionality is broken

As we add more assertions and specific test cases later, we'll use these more descriptive names to make it clear what each test is verifying.

### 1.2. Run the test

Let's run the test to see what happens.

```bash
nf-test test tests/main.nf.test
```

```console title="nf-test pipeline fail"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures' FAILED (4.652s)

  Assertion failed:

  assert workflow.success
         |        |
         workflow false

  Nextflow stdout:

  ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv

   -- Check '/workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/meta/nextflow.log' file for details
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.679s (1 failed)
```

The test fails! What happened?

1. nf-test tried to run the pipeline as is, using the settings in the `when` block:

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

2. nf-test checked the status of the pipeline and compared it to the `when` block:

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

Note how nf-test has reported the pipeline failed and provided the error message from Nextflow:

```console title="Error"
ERROR ~ No such file or directory: /workspaces/training/side-quests/nf-test/.nf-test/tests/693ba951a20fec36a5a9292ed1cc8a9f/greetings.csv
```

So what was the issue? Remember the pipeline has a greetings.csv file in the project directory. When nf-test runs the pipeline, it will look for this file, but it can't find it. The file is there, what's happening? Well, if we look at the path we can see the test is occurring in the path `./nf-test/tests/longHashString/`. Just like Nextflow, nf-test creates a new directory for each test to keep everything isolated. The data file is not located in there so we must correct the path to the file in the original test.

Let's go back to the test file and change the path to the file in the `when` block.

You may be wondering how we're going to point to the root of the pipeline in the test. Since this is a common situation, nf-test has a range of global variables that we can use to make our lives easier. You can find the full list [here](https://www.nf-test.com/docs/testcases/global_variables/) but in the meantime we'll use the `projectDir` variable, which means the root of the pipeline project.

_Before:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3 4"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

_After:_

```groovy title="tests/main.nf.test" linenums="1" hl_lines="3"
when {
    params {
        input_file = "${projectDir}/greetings.csv"
    }
}
```

Let's run the test again to see if it works.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.619s)


SUCCESS: Executed 1 tests in 1.626s
```

Success! The pipeline runs successfully and the test passes. Run it as many times as you like and you will always get the same result!

By default, the Nextflow output is hidden, but to convince yourself that nf-test is definitely running the workflow, you can use the `--verbose` flag:

```bash
nf-test test tests/main.nf.test --verbose
```

```console title="Pipeline runs all processes"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [693ba951] 'Should run without failures'
    > Nextflow 24.10.4 is available - Please consider updating your version to it
    > N E X T F L O W  ~  version 24.10.0
    > Launching `/workspaces/training/side-quests/nf-test/main.nf` [zen_ampere] DSL2 - revision: bbf79d5c31
    > [2b/61e453] Submitted process > sayHello (2)
    > [31/4e1606] Submitted process > sayHello (1)
    > [bb/5209ee] Submitted process > sayHello (3)
    > [83/83db6f] Submitted process > convertToUpper (2)
    > [9b/3428b1] Submitted process > convertToUpper (1)
    > [ca/0ba51b] Submitted process > convertToUpper (3)
    PASSED (5.206s)


SUCCESS: Executed 1 tests in 5.239s
```

### 1.3. Add assertions

A simple check is to ensure our pipeline is running all the processes we expect and not skipping any silently. Remember our pipeline runs 6 processes, one called `sayHello` and one called `convertToUpper` for each of the 3 greetings.

Let's add an assertion to our test to check the pipeline runs the expected number of processes. We'll also update our test name to better reflect what we're testing.

**Before:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1"
    test("Should run without failures") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
        }

    }
```

**After:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="1 11"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

The test name now better reflects what we're actually verifying - not just that the pipeline runs without failing, but that it runs the expected number of processes.

Let's run the test again to see if it works.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with assertions"
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run successfully with correct number of processes' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

Success! The pipeline runs successfully and the test passes. Now we have began to test the details of the pipeline, as well as the overall status.

### 1.4. Test the output

Let's add an assertion to our test to check the output file was created. We'll add it as a separate test, with an informative name, to make the results easier to interpret.

**Before:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }
```

**After:**

```groovy title="tests/main.nf.test" linenums="1" hl_lines="14-33"
    test("Should run successfully with correct number of processes") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert workflow.success
            assert workflow.trace.tasks().size() == 6
        }

    }

    test("Should produce correct output files") {

        when {
            params {
                input_file = "${projectDir}/greetings.csv"
            }
        }

        then {
            assert file("$launchDir/results/Bonjour-output.txt").exists()
            assert file("$launchDir/results/Hello-output.txt").exists()
            assert file("$launchDir/results/Holà-output.txt").exists()
            assert file("$launchDir/results/UPPER-Bonjour-output.txt").exists()
            assert file("$launchDir/results/UPPER-Hello-output.txt").exists()
            assert file("$launchDir/results/UPPER-Holà-output.txt").exists()
        }

    }
```

Run the test again to see if it works.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with file assertions"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)


SUCCESS: Executed 2 tests in 15.165s
```

Success! The tests pass because the pipeline completed successfully, the correct number of processes ran and the output files were created. This should also show you how useful it is to provide those informative names for your tests.

This is just the surface, we can keep writing assertions to check the details of the pipeline, but for now let's move on to testing the internals of the pipeline.

### Takeaway

You know how to write an nf-test for a pipeline.

### What's next?

Learn how to test a Nextflow process.

---

## 2. Test a Nextflow process

We don't have to write tests for every part of the pipeline, but the more tests we have the more comprehensive we can be about the pipeline and the more confident we can be that it's working as expected. In this section we're going to test both processes in the pipeline as individual units.

### 2.1. Test the `sayHello` process

Let's start with the `sayHello` process.

Let's use the `nf-test generate` command again to generate tests for the process.

```bash
nf-test generate process main.nf
```

```console title="Output"
> nf-test generate process main.nf

Load source file '/workspaces/training/side-quests/nf-test/main.nf'
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.sayhello.nf.test
Wrote process test file '/workspaces/training/side-quests/nf-test/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

Let's focus for now on the `sayhello` process in the `main.sayhello.nf.test` file.

Let's open the file and take a look at the contents.

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {

    name "Test Process sayHello"
    script "main.nf"
    process "sayHello"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

As before, we start with the test details, followed by the `when` and `then` blocks. However, we also have an additional `process` block which allows us to define the inputs to the process.

Let's run the test to see if it works.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [1eaad118] 'Should run without failures' FAILED (4.876s)

  Assertion failed:

  assert process.success
         |       |
         |       false
         sayHello

  Nextflow stdout:

  Process `sayHello` declares 1 input channel but 0 were specified

   -- Check script '/workspaces/training/side-quests/nf-test/.nf-test-1eaad118145a1fd798cb07e7dd75d087.nf' at line: 38 or see '/workspaces/training/side-quests/nf-test/.nf-test/tests/1eaad118145a1fd798cb07e7dd75d087/meta/nextflow.log' file for more details
  Nextflow stderr:

FAILURE: Executed 1 tests in 4.884s (1 failed)
```

The test fails because the `sayHello` process declares 1 input channel but 0 were specified. Let's fix that by adding an input to the process. Remember from [Hello Workflow](../hello_nextflow/03_hello_workflow.md) (and the warmup section above) that our `sayHello` process takes a single value input, which we will need to provide. We should also fix the test name to better reflect what we're testing.

**Before:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**After:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Let's run the test again to see if it works.

```console title="nf-test pipeline pass"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.604s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.611s
```

Success! The test passes because the `sayHello` process ran successfully and the output was created.

### 2.2. Check out the snapshot created by the test

If we look at the `tests/main.sayhello.nf.test` file, we can see it uses a method `snapshot()` in the assertion block:

```groovy title="tests/main.sayhello.nf.test"
assert snapshot(process.out).match()
```

This is telling nf-test to create a snapshot of the output of the `sayHello` process. Let's take a look at the contents of the snapshot file.

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

We won't print it here, but you should see a JSON file containing details of the process and process outputs. In particular, we can see a line that looks like this:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

This represents the outputs created by the `sayHello` process, which we are testing explicitly. If we re-run the test, the program will check that the new output matches the output that was originally recorded. This is a quick, simple way of testing that process outputs don't change, which is why nf-test provides it as a default.

!!!warning

    That means we have to be sure that the output we record in the original run is correct!

If, in the course of future development, something in the code changes that causes the output to be different, the test will fail and we will have to determine whether the change is expected or not.

- If it turns out that something in the code broke, we will have to fix it, with the expectation that the fixed code will pass the test.
- If it is an expected change (e.g., the tool has been improved and the results are better) then we will need to update the snapshot to accept the new output as the reference to match. nf-test has a parameter `--update-snapshot` for this purpose.

We can run the test again and see the test should pass:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and produce correct output' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Success! The test passes because the `sayHello` process ran successfully and the output matched the snapshot.

### 2.3. Alternative to Snapshots: Direct Content Assertions

While snapshots are great for catching any changes in output, sometimes you want to verify specific content without being so strict about the entire file matching. For example:

- When parts of the output might change (timestamps, random IDs, etc.) but certain key content must be present
- When you want to check for specific patterns or values in the output
- When you want to make the test more explicit about what constitutes success

Here's how we could modify our test to check specific content:

**Before:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 6 17"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**After:**

```groovy title="tests/main.sayhello.nf.test" linenums="1" hl_lines="1 5 16 17"
    test("Should run without failures and contain expected greeting") {

        when {
            params {
                // define parameters here
            }
            process {
                """
                input[0] = "hello"
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('hello')
            assert !path(process.out[0][0]).readLines().contains('HELLO')
        }

    }
```

Note that nf-test sees the process outputs as a list of lists, so `process.out[0][0]` is fetching the first part of the first channel item (or 'emission') from this process.

This approach:

- Makes it clear exactly what we expect in the output
- Is more resilient to irrelevant changes in the output
- Provides better error messages when tests fail
- Allows for more complex validations (regex patterns, numerical comparisons, etc.)

Let's run the test to see if it works.

```bash title="nf-test pipeline pass"
nf-test test tests/main.sayhello.nf.test
```

```console title="Process test fails"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [58df4e4b] 'Should run without failures and contain expected greeting' PASSED (7.196s)


SUCCESS: Executed 1 tests in 7.208s
```

### 2.4. Test the `convertToUpper` process

Let's open the `tests/main.converttoupper.nf.test` file and take a look at the contents:

```groovy title="tests/main.converttoupper.nf.test"
nextflow_process {

    name "Test Process convertToUpper"
    script "main.nf"
    process "convertToUpper"

    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

}
```

This is a similar test to the `sayHello` process, but it's testing the `convertToUpper` process. We know this one will fail because just like with `sayHello`, the `convertToUpper` process takes a single path input, but we haven't specified one.

We now need to supply a single input file to the convertToUpper process, which includes some text that we want to convert to uppercase. There are lots of ways we could do this:

- We could create a dedicated file to test
- We could re-use the existing data/greetings.csv file
- We could create it on the fly within the test

For now, let's re-use the existing data/greetings.csv file using the example we used with the pipeline level test. As before, we can name the test to better reflect what we're testing, but this time let's leave it to 'snapshot' the content rather than checking for specific strings (as we did in the other process).

**Before:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10 11"
    test("Should run without failures") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                // define inputs of the process here. Example:
                // input[0] = file("test-file.txt")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

**After:**

```groovy title="tests/main.converttoupper.nf.test" linenums="1" hl_lines="1 10"
    test("Should run without failures and produce correct output") {

        when {
            params {
                // define parameters here. Example:
                // outdir = "tests/results"
            }
            process {
                """
                input[0] = "${projectDir}/greetings.csv"
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

And run the test!

```bash title="nf-test pipeline pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.755s)
  Snapshots:
    1 created [Should run without failures and produce correct output]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 1.764s
```

Note, we have created a snapshot file for the `convertToUpper` process at `tests/main.converttoupper.nf.test.snap`. If we run the test again, we should see the nf-test passes again.

```bash title="nf-test process convertToUpper pass"
nf-test test tests/main.converttoupper.nf.test
```

```console title="nf-test process convertToUpper pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### Takeaway

You know how to write tests for a Nextflow process and run them.

### What's next?

Learn how to run tests for everything at once!

## 3. Run tests for the entire repository

Running nf-test on each component is fine, but laborious and error prone. Can't we just test everything at once?

Yes we can!

Let's run nf-test on the entire repo.

### 3.1. Run nf-test on the entire repo

We can run nf-test on the entire repo by running the `nf-test test` command.

```bash
nf-test test .
```

Note, we are just using the `.` to run everything from our current directory. This will include every test!

```console title="nf-test repo pass"
> nf-test test tests/main.converttoupper.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures and produce correct output' PASSED (1.798s)

Test Workflow main.nf

  Test [f0e08a68] 'Should run successfully with correct number of processes' PASSED (8.144s)
  Test [d7e32a32] 'Should produce correct output files' PASSED (6.994s)

Test Process sayHello

  Test [f91a1bcd] 'Should run without failures and contain expected greeting' PASSED (1.664s)


SUCCESS: Executed 3 tests in 5.007s
```

Check that out! We ran 3 tests, 1 for each process and 1 for the whole pipeline with a single command. Imagine how powerful this is on a large codebase!

---

## Summary

In this side quest, you've learned to leverage nf-test's features to create and run tests for individual processes as well as end-to-end tests for the entire pipeline.
You're now aware of the main two approaches to output validation, snapshots and direct content assertions, and and when to use either one.
You also know how to run tests either one by one or for an entire project.

Applying these techniques in your own work will enable you to ensure that:

- Your code works as expected
- Changes don't break existing functionality
- Other developers can contribute with confidence
- Problems can be identified and fixed quickly
- Output content matches expectations

### Key patterns

<!-- TODO: Can we add snippets of code below to illustrate? -->

1. Pipeline-level tests:
   - Basic success testing
   - Process count verification
   - Output file existence checks
2. Process-level tests
3. Two approaches to output validation:
   - Using snapshots for complete output verification
   - Using direct content assertions for specific content checks
4. Running all tests in a repository with a single command

### Additional resources

Check out the [nf-test documentation](https://www.nf-test.com/) for more advanced testing features and best practices. You might want to:

- Add more comprehensive assertions to your tests
- Write tests for edge cases and error conditions
- Set up continuous integration to run tests automatically
- Learn about other types of tests like workflow and module tests
- Explore more advanced content validation techniques

**Remember:** Tests are living documentation of how your code should behave. The more tests you write, and the more specific your assertions are, the more confident you can be in your pipeline's reliability.

---

## What's next?

Return to the [menu of Side Quests](./index.md) or click the button in the bottom right of the page to move on to the next topic in the list.
