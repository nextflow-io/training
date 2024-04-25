# Part 4: Hello nf-test

[short blurb about testing]

For a more detailed overview, read the [blog post about nf-test](https://nextflow.io/blog/2024/nf-test-in-nf-core.html) on the nf-core blog.

---

## 0. Warmup

We start from a base workflow called `hello-nf-test.nf`, which corresponds to the workflow we produced in Part 3: Hello modules (equivalent to `scripts/hello-modules-3.nf`). 

This is a modularized pipeline; the processes are in local modules and the parameter declarations are in a configuration file. If you completed the previous parts of the training course, then you already have everything you need in the working directory. However, if you're picking this up here, you need to copy the `nextflow.config` file and the `modules` folder from `scripts/` to the working directory.

### 0.1 Run the workflow to verify that it produces the expected outputs 

```bash
nextflow run hello-nf-test.nf
```

[TODO: as before, briefly describe outputs; refer to Part 2 for full details]

### 0.2 Initialize nf-test

Run the following command in the terminal:

```bash
nf-test init
```

This should produce the following output:

```bash
🚀 nf-test 0.8.4
https://code.askimed.com/nf-test
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Project configured. Configuration is stored in nf-test.config
```

It also creates a `tests` directory containing a configuration file stub.

---

## 1. Test a process for success and matching outputs

We're going to start by adding a test for the `SAMTOOLS_INDEX` process. It's a very simple process that takes a single input file (which we have test data for on hand) and generates an index file. We want to test that the process runs successfully and that the file it produces is always the same for a given input.

### 1.1 Run the `nf-test generate` command to create a test file stub

```bash
nf-test generate process modules/local/samtools/index/main.nf
```

This creates a test file stub created in same directory as `main.nf`, summarized in the terminal output as follows:

```bash
Load source file '/workspace/gitpod/hello-nextflow/modules/local/samtools/index/main.nf'
Wrote process test file '/workspace/gitpod/hello-nextflow/tests/modules/local/samtools/index/main.nf.test

SUCCESS: Generated 1 test files.
```

You can navigate to the directory in the file explorer and open the file, which should contain the following code:

```groovy
nextflow_process {

    name "Test Process SAMTOOLS_INDEX"
    script "modules/local/samtools/index/main.nf"
    process "SAMTOOLS_INDEX"

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

[TODO: Add a brief explanation of the two `assert` statements]

As you can see, this is just a stub, it's not yet a functional test. We need to add parameters and inputs that will be fed to the process.

### 1.2 Move the test file and update the script path

Before we get to work on filling out the test, we need to move the file to its definitive location. The preferred convention is to ship tests in a `tests` directory co-located with each module's `main.nf` file, so we create a `tests/` directory in the module's directory:

```bash
mkdir -p modules/local/samtools/index/tests
```

Then we move the test file there:

```bash
mv tests/modules/local/samtools/index/main.nf.test modules/local/samtools/index/tests/
```

As a result, we can update the full path in the `script` section of the test file to a relative path:

_Before:_
```groovy
name "Test Process SAMTOOLS_INDEX"
script "modules/local/samtools/index/main.nf"
process "SAMTOOLS_INDEX"
```

_After:_
```groovy
name "Test Process SAMTOOLS_INDEX"
script "../main.nf"
process "SAMTOOLS_INDEX"
```


### 1.3 Provide inputs to the test process

The stub file includes a placeholder that we need to replace with an actual test input:

_Before:_
```groovy
process {
    """
    // define inputs of the process here. Example:
    // input[0] = file("test-file.txt")
    """
}
```

We choose one of our available data files and plug it into the `process` block:

_After:_
```groovy
process {
    """
    input[0] = [ [id: 'NA12882' ], file("/workspace/gitpod/hello-nextflow/data/bam/reads_son.bam") ]
    """
}
```

### 1.4 Rename the test based on the primary test input

The stub file gives the test a generic name referring to the assertion that it should run without failures. Since we added a specific input, it's good practice to rename the test accordingly.

_Before:_
```groovy
test("Should run without failures") {
```

This takes an arbitrary string, so we could put anything we want. Here we choose to refer to the file name and its format:

_After:_
```groovy
test("reads_son [bam]") {
```

### 1.5 Specify test parameters

The `params` block in the stub file includes a placeholder for parameters:

_Before:_
```groovy
params {
    // define parameters here. Example:
    // outdir = "tests/results"
}
```

We use it to specify a location for the results to be output, using the default suggestion:

_After:_
```groovy
params {
    outdir = "tests/results"
}
```

### 1.6 Run the test and examine the output

```bash
nf-test test modules/local/samtools/index/tests/main.nf.test
```

This should produce the following output: 

```bash
🚀 nf-test 0.8.4
https://code.askimed.com/nf-test
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process SAMTOOLS_INDEX

  Test [bc664c47] 'reads_son [bam]' PASSED (10.06s)
  Snapshots:
    1 created [reads_son [bam]]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 10.068s
```

The test verified the first assertion, that the process should complete successfully.

Additionally, this also produces a snapshot file called `main.nf.test.snap` that captures all the output channels and the MD5SUMs of all elements. If we re-run the test, the program will check that the new output matches the output that was originally recorded. 

Note: That does mean that we have to be sure that the output we record in the original run is correct.

If, in the course of future development, something in the code changes that causes the output to be different, the test will fail and we will have to determine whether the change is expected or not. If it turns out that something in the code broke, we will have to fix it, with the expectation that the fixed code will pass the test. If it is an expected change (e.g. the tool has been improved and the results are better) then we will need to update the snapshot to accept the new output as the reference to match, using the parameter `--update-snapshot` when we run the test command. 

### 1.7 Add more tests to `SAMTOOLS_INDEX`

Sometimes it's useful to test a range of different input files to ensure we're testing for a variety of potential issues. We can add as many tests as we want inside the same test file for a module.

Try adding the following tests to the module's test file:

```groovy
test("reads_mother [bam]") {

    when {
        params {
            outdir = "tests/results"
        }
        process {
            """
            input[0] = [ [id: 'NA12878' ], file("/workspace/gitpod/hello-nextflow/data/bam/reads_mother.bam") ]
            """
        }
    }

    then {
        assert process.success
        assert snapshot(process.out).match()
    }

}
```

And:

```groovy
test("reads_father [bam]") {

    when {
        params {
            outdir = "tests/results"
        }
        process {
            """
            input[0] = [ [id: 'NA12877' ], file("/workspace/gitpod/hello-nextflow/data/bam/reads_father.bam") ]
            """
        }
    }

    then {
        assert process.success
        assert snapshot(process.out).match()
    }

}
```

These simply go one after another in the test file.

Note: Watch those curly braces, make sure they're all paired up appropriately...

### 1.8 Run the test suite and update the snapshot

```bash
nf-test test modules/local/samtools/index/tests/main.nf.test --update-snapshot
```

This should produce the following output:

```bash
🚀 nf-test 0.8.4
https://code.askimed.com/nf-test
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Warning: every snapshot that fails during this test run is re-record.

Test Process SAMTOOLS_INDEX

  Test [bc664c47] 'reads_son [bam]' PASSED (9.6s)
  Test [f413ec92] 'reads_mother [bam]' PASSED (9.138s)
  Test [99a73481] 'reads_father [bam]' PASSED (9.536s)


SUCCESS: Executed 3 tests in 28.281s
```

Note the warning, referring to the effect of the `--update-snapshot` parameter.

[Anything to add? Maybe a note about being careful to keep unit tests light by using small data, otherwise runtime adds up. Here we could be using even smaller data snippets if we cared about reducing runtime for a production environment.]

---

## 2. Add tests to a chained process and test for contents

Now that we know how to handle the simplest case, we're going to kick things up a notch with the `GATK_HAPLOTYPECALLER` process. As the second step in our pipeline, its input depends on the output of another process. We can deal with this in two ways: either manually generate some static test data that is suitable as intermediate input to the process, or we can use a special [setup method](https://www.nf-test.com/docs/testcases/setup/) to handle it dynamically for us.

Spoiler: we're going to use the setup method.

We're going to test that the process runs successfully and we're also going to check for specific content within the output file. 

### 2.1 Generate the test file stub and move it to the right place

As previously, first we generate the file stub:

```bash
nf-test generate process modules/local/gatk/haplotypecaller/main.nf
```

Then we create a directory for it co-located with the module's `main.nf` file:

```bash
mkdir -p modules/local/gatk/haplotypecaller/tests
```

And move the test stub file there:

```bash
mv tests/modules/local/gatk/haplotypecaller/main.nf.test modules/local/gatk/haplotypecaller/tests/
```

Finally, don't forget to update the script path:

_Before:_
```groovy
name "Test Process GATK_HAPLOTYPECALLER"
script "modules/local/gatk/haplotypecaller/main.nf"
process "GATK_HAPLOTYPECALLER"
```

_After:_
```groovy
name "Test Process GATK_HAPLOTYPECALLER"
script "../main.nf"
process "GATK_HAPLOTYPECALLER"
```

### 2.2 Provide inputs using the setup method 




### 3.3 [run test and examine output]

### 3.4 [run again and observe failure]

### 3.5 [use contains assertion method]

### 3.6 [run and observe success]

---

## 4. Add more inputs to second module test

[...]

### 4.x [run full set of tests]

---

## 5. Add local input tests to third process

[ensure test files are in dir]

### 5.1 [generate process test files for third module]

### 5.2 [set up test (may be more steps?)]

### 5.3 [run test]

---

## 6. Generate pipeline level test

### 6.1 [generate pipeline-level nf-test file]

### 6.2 [run test (assert success)]

### 6.3 [optionally add more specific assertions if needed]





