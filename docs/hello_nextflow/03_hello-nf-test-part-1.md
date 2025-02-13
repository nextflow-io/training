# Part 7: Hello nf-test

Being able to systematically test that every part of your workflow is doing what it's supposed to do is critical for reproducibility and long-term maintenance.
And it's also helpful during the development process!

This is (hopefully) not controversial, but let's take a minute to talk about why testing is so important.

Let's imagine you're developing a workflow. One of the first things you will do is grab some test data that you know is valid and should produce a result. You add the first process to the pipeline and wire it up to your inputs to make it work. Then, to check it's all working, you run it on the test data. Assuming it works, you move on to the next process and run the test data again. You repeat this process until you have a pipeline that you're happy with.

Maybe you add a simple true or false parameter such as `--skip_process`, now you must run the pipeline twice, once with each parameter to make sure it works as expected. But wait, how do we check if the `--skip_process` actually skips the process? We have to dig into the outputs or check the log files! This is a pain and prone to error.

As you develop your pipeline, it will quickly become so complex that manually testing every iteration is slow and error prone. Furthermore, if you do find an error it will be very difficult to pin down exactly where in your pipeline the error is coming from. This is where unit testing comes in.

Testing allows you to systematically check that every part of your pipeline is working as expected. The benefits to a developer are huge:

-   **Confidence**: Because the tests cover the entire pipeline, you can be be confident changing something doesn't affect anything else
-   **Trust**: When multiple developers work on the pipeline, they know the other developers haven't broken the pipeline and every component.
-   **Transparency**: The tests show where a pipeline is failing and make it easier to track down the problem. They also function as a form of documentation, showing how to run a process or workflow.
-   **Speed**: Because the tests are automated, they can be run very quickly and repeatedly. You can iterate quickly with less fear of introducing new bugs.

There are lots of different types of tests we can write:

1. **Module-level tests**: For individual processes
2. **Workflow-level tests**: For a single workflow
3. **Pipeline-level tests**: For the pipeline as a whole
4. **Integration tests**: For the interaction of the pipeline with other systems
5. **Performance tests**: For the speed and efficiency of the pipeline
6. **Stress tests**: To identify the limits of the pipeline

[**nf-test**](https://www.nf-test.com/) is a tool that allows you to write module, workflow and pipeline level test. In short, it allows you to systematically check every individual part of the pipeline is working as expected, _in isolation_.

In this part of the training, we're going to show you how to use nf-test to write module-level tests for the three processes in our pipeline.

---

## 0. Warmup

Let's move into the project directory.
If you're continuing on directly from Part 2, you'll need to move up one directory first.

```bash
cd hello-nf-test-part1
```

The `hello-nf-test-part1` directory has the same content and structure that you're expected to end up with in `hello-world` on completion of Part 1.

```console title="Directory contents"
hello-nf-test-part1/
├── data
│   └── greetings.csv
└── main.nf
```

For a detailed description of the files, see the warmup section in Part 1.

### 0.1. Run the workflow

Let's run the workflow to make sure it's working as expected.

```bash
nextflow run main.nf
```

CONGRATULATIONS! You just ran a unit test!

"Wait, what? I just ran the workflow and it worked! How is that a unit test?"

Good question!

Let's break down what just happened.

You ran the workflow with the default parameters, you confirmed it worked and you're happy with the results. This is the essence of unit testing. As you go through this course you will notice we always run the workflow at the start to confirm everything is set up correctly. This gives us the confidence to continue with the training, knowing the code started in a working state.

Unit testing software essentially does this process for us. Let's replace our simple `nextflow run main.nf` with a standardised test provided by nf-test.

### Takeaway

You should be able to 'test' a pipeline by manually running it.

### What's next?

Initialize `nf-test`.

---

## 1.0. Initialize `nf-test`

The `nf-test` package provides an initialization command that sets up a few things in order for us to start developing tests for our project.

```bash
nf-test init
```

This should produce the following output:

```bash
🚀 nf-test 0.9.1
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

It also creates a `tests` directory containing a configuration file stub.

## 1.1. Generate an nf-test

`nf-test` comes with a set of tools for building nf-test files, saving us the majority of the work. These come under the subcommand `generate`. Let's generate a test for the pipeline:

```bash
nf-test generate pipeline main.nf
```

```console title="Output"
> nf-test generate pipeline main.nf

Load source file '/workspace/gitpod/hello-nextflow/hello-nf-test-part1/main.nf'
Wrote pipeline test file '/workspace/gitpod/hello-nextflow/hello-nf-test-part1/tests/main.nf.test

SUCCESS: Generated 1 test files.
```

This will create a `tests` directory with a `main.nf.test` and a `nextflow.config` file. This is our pipeline level test file. If you run `tree tests/` you should see something like this:

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

-   `name`: The name of the test.
-   `script`: The path to the pipeline script.

The `test` block is the actual test. It contains the following:

-   `when`: The conditions under which the test should be run. This includes the parameters that will be used to run the pipeline.
-   `then`: The assertions that should be made. This includes the expected outcomes of the pipeline.

In plain English, the logic of the test reads as follows:
"**When** these _parameters_ are provided to this _pipeline_, **then** we expect to see these results."

This isn't a functional test, we will demonstrate how to turn it into one in the next section.

## 1.2. Run the test

Let's run the test to see what happens.

```console title="nf-test pipeline fail"
> nf-test test tests/main.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' FAILED (1.325s)

  Assertion failed:

  assert workflow.success
         |        |
         workflow false

  Nextflow stdout:

  ERROR ~ No such file or directory: /workspace/gitpod/hello-nextflow/hello-nf-test-part1/.nf-test/tests/1d4aaf120a1c8314a4310f6a54a0c243/data/greetings.csv

   -- Check '/workspace/gitpod/hello-nextflow/hello-nf-test-part1/.nf-test/tests/1d4aaf120a1c8314a4310f6a54a0c243/meta/nextflow.log' file for details
  Nextflow stderr:
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
ERROR ~ No such file or directory: /workspace/gitpod/hello-nextflow/hello-nf-test-part1/.nf-test/tests/1d4aaf120a1c8314a4310f6a54a0c243/data/greetings.csv
```

So what was the issue? Remember the pipeline has a greetings.csv file in the `data` directory. When nf-test runs the pipeline, it will look for this file, but it can't find it. The file is there, what's happening? Well, if we look at the path we can see the test is occurring in the path `./nf-test/tests/longHashString/`. Just like Nextflow, nf-test creates a new directory for each test to keep everything isolated. The data file is not located in there so we must correct the path to the file in the original test.

Let's go back to the test file and change the path to the file in the `when` block.

You may be wondering how we're going to point to the root of the pipeline in the test. Since this is a common situation, nf-test has a range of global variables that we can use to make our lives easier. You can find the full list [here](https://www.nf-test.com/docs/testcases/global_variables/) but in the meantime we'll use the `projectDir` variable, which means the root of the pipeline project.

_Before:_

```groovy title="tests/main.nf.test"
when {
    params {
        // define parameters here. Example:
        // outdir = "tests/results"
    }
}
```

_After:_

```groovy title="tests/main.nf.test"
when {
    params {
        input_file = "${projectDir}/data/greetings.csv"
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

## 1.3. Add assertions

A simple check is to ensure our pipeline is running all the processes we expect and not skipping any silently. Remember our pipeline runs 6 processes, one called `sayHello` and one called `convertToUpper` for each of the 3 greetings.

Let's add an assertion to our test to check the pipeline runs the expected number of processes.

**Before:**

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
}
```

**After:**

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
    assert workflow.trace.tasks().size() == 6
}
```

The `workflow.trace` object includes information about the pipeline which we can check. In this case, we're checking the number of tasks is correct.

Let's run the test again to see if it works.

```bash title="nf-test pipeline pass"
nf-test test tests/main.nf.test
```

```console title="Pipeline passes with assertions"
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.567s)


SUCCESS: Executed 1 tests in 1.588s
```

Success! The pipeline runs successfully and the test passes. Now we have began to test the details of the pipeline. as well as the overall status.

## 1.4. Test the output

Let's add an assertion to our test to check the output file was created.

**Before:**

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
    assert workflow.trace.tasks().size() == 6
}
```

**After:**

```groovy title="tests/main.nf.test"
then {
    assert workflow.success
    assert workflow.trace.tasks().size() == 6
    assert file("results/Bonjour-output.txt").exists()
    assert file("results/Hello-output.txt").exists()
    assert file("results/Holà-output.txt").exists()
    assert file("results/UPPER-Bonjour-output.txt").exists()
    assert file("results/UPPER-Hello-output.txt").exists()
    assert file("results/UPPER-Holà-output.txt").exists()
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

  Test [1d4aaf12] 'Should run without failures' PASSED (1.591s)


SUCCESS: Executed 1 tests in 1.612s
```

Success! The test passes because the pipeline completed successfully, the correct number of processes ran and the output files were created.

This is just the surface, we can keep writing assertions to check the details of the pipeline, but for now let's move on to testing the internals of the pipeline.

### Takeaway

You know how to write an nf-test for a pipeline.

### What's next?

Learn how to test a Nextflow process.

---

## 2 Test a Nextflow process

We don't have to write tests for every part of the pipeline, but the more tests we have the more comprehensive we can be about the pipeline and the more confident we can be that it's working as expected. In this section we're going to test both processes in the pipeline.

### 2.1. Test the `sayHello` process

Let's start with the `sayHello` process.

Let's use the `nf-test generate` command again to generate tests for the entire pipeline.

```bash
nf-test generate process main.nf
```

```console title="Output"
> nf-test generate process main.nf

Load source file '/workspace/gitpod/hello-nextflow/hello-nf-test-part1/main.nf'
Wrote process test file '/workspace/gitpod/hello-nextflow/hello-nf-test-part1/tests/main.sayhello.nf.test
Wrote process test file '/workspace/gitpod/hello-nextflow/hello-nf-test-part1/tests/main.converttoupper.nf.test

SUCCESS: Generated 2 test files.
```

Note we have created a test for the `converttoupper` process as well. We can ignore that for now and focus on the `sayhello` process in the `main.sayhello.nf.test` file.

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

  Test [f91a1bcd] 'Should run without failures' FAILED (1.347s)

  Assertion failed:

  assert process.success
         |       |
         |       false
         sayHello

  Nextflow stdout:

  Process `sayHello` declares 1 input channel but 0 were specified

   -- Check script '/workspace/gitpod/hello-nextflow/hello-nf-test-part1/.nf-test-f91a1bcdafbf4c8bbd32858acdd8afd2.nf' at line: 38 or see '/workspace/gitpod/hello-nextflow/hello-nf-test-part1/.nf-test/tests/f91a1bcdafbf4c8bbd32858acdd8afd2/meta/nextflow.log' file for more details
  Nextflow stderr:

FAILURE: Executed 1 tests in 1.35s (1 failed)
```

The test fails because the `sayHello` process declares 1 input channel but 0 were specified. Let's fix that by adding an input to the process. Remember from part 1, our `sayHello` process takes a single value input.

**Before:**

```groovy title="tests/main.sayhello.nf.test"
process {
    """
    // define inputs of the process here. Example:
    // input[0] = file("test-file.txt")
    """
}
```

**After:**

```groovy title="tests/main.sayhello.nf.test"
process {
    """
    input[0] = "hello"
    """
}
```

Let's run the test again to see if it works.

```console title="nf-test pipeline pass"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures' PASSED (1.604s)
  Snapshots:
    1 created [Should run without failures]


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

This is a snapshot of the output of the `sayHello` process. Let's take a look at the contents of the snapshot file.

```console title="Snapshot file contents"
code tests/main.sayhello.nf.test.snap
```

We won't print it here, but you should see a JSON file containing details of the process and process outputs. In particular, we can see a line that looks like this:

```json title="Snapshot file contents"
"0": [
    "hello-output.txt:md5,b1946ac92492d2347c6235b4d2611184"
]
```

This represents the outputs created by the `sayHello` process, which we are testing explicitly. If we re-run the test, the program will check that the new output matches the output that was originally recorded. This is a quick, simple way of testing process outputs which is why nf-test provides it as a default.

!!!warning

    That means we have to be sure that the output we record in the original run is correct!

If, in the course of future development, something in the code changes that causes the output to be different, the test will fail and we will have to determine whether the change is expected or not.

-   If it turns out that something in the code broke, we will have to fix it, with the expectation that the fixed code will pass the test.
-   If it is an expected change (e.g., the tool has been improved and the results are better) then we will need to update the snapshot to accept the new output as the reference to match, using the parameter `--update-snapshot` when we run the test command.

For now though, we can run the test again and see the test should pass:

```console title="nf-test process pass with snapshot"
> nf-test test tests/main.sayhello.nf.test

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process sayHello

  Test [f91a1bcd] 'Should run without failures' PASSED (1.675s)


SUCCESS: Executed 1 tests in 1.685s
```

Success! The test passes because the `sayHello` process ran successfully and the output matched the snapshot.

### 2.3. Test the `convertToUpper` process

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

-   We could create a dedicated file to test
-   We could re-use the existing data/greetings.csv file
-   We could create it on the fly within the test

For now, let's re-use the existing data/greetings.csv file using the example we used with the pipeline level test.

**Before:**

```groovy title="tests/main.converttoupper.nf.test"
process {
    """
    // define inputs of the process here. Example:
    // input[0] = file("test-file.txt")
    """
}
```

**After:**

```groovy title="tests/main.converttoupper.nf.test"
process {
    """
    input[0] = "${projectDir}/data/greetings.csv"
    """
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

  Test [c59b6044] 'Should run without failures' PASSED (1.755s)
  Snapshots:
    1 created [Should run without failures]


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

  Test [c59b6044] 'Should run without failures' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s
```

### Takeaway

You know how to write tests for a Nextflow process and run them.

### What's next?

Learn how to run tests for everything at once!

## 3. Run tests for the entire repo

Running nf-test on each component is fine, but laborious and error prone. Can't we just test a test on everything at once?

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

  Test [c59b6044] 'Should run without failures' PASSED (1.798s)


SUCCESS: Executed 1 tests in 1.811s

> nf-test test .

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process convertToUpper

  Test [c59b6044] 'Should run without failures' PASSED (1.663s)

Test Workflow main.nf

  Test [1d4aaf12] 'Should run without failures' PASSED (1.652s)

Test Process sayHello

  Test [f91a1bcd] 'Should run without failures' PASSED (1.664s)


SUCCESS: Executed 3 tests in 5.007s
```

Check that out! We ran 3 tests, 1 for each process and 1 for the whole pipeline with a single command. Imagine how powerful this is on a large codebase!

### Takeaway

You know how to run tests for the entire repo with a single command.

### What's next?

This is a lot to learn, so celebrate and take a big break!

If you already know about Nextflow and wish to learn more about nf-test, skip to nf-test part 2 which uses more real world examples.

If you are still at the beginning of your Nextflow journey, you can continue to learn about containers.
