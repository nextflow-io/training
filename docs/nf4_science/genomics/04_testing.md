# Part 4: Adding tests

In the first part of this course, you built a variant calling pipeline that was completely linear and processed each sample's data independently of the others.

In the second part, we showed you how to use channels and channel operators to implement joint variant calling with GATK.

In the third part, we modularized the pipeline.

In this part of the training, we're going to show you how to use [**nf-test**](https://www.nf-test.com/), a testing framework that integrates well with Nextflow and makes it straightforward to add both module-level and workflow-level tests to your pipeline. To follow this part of the training, you should have completed Part 1, Part 2, and Part 3, as well as the [nf-test side quest](../../side_quests/nf-test.md), which covers the basics of nf-test, and why testing is important.

---

## 0. Warmup

!!! note

     Make sure you're in the correct working directory:
     `cd /workspaces/training/nf4-science/genomics`

If you worked through the previous parts of this training course, you should have a working version of the genomics pipeline, with a modules directory structure like:

```console title="Directory structure"
modules/
├── gatk
│   ├── haplotypecaller
│   │   └── main.nf
│   └── jointgenotyping
│       └── main.nf
└── samtools
    └── index
        └── main.nf
```

This modules directory can be found in the `solutions` directory if you need it.

We're going to start with the same workflow as in Part 3, which we've provided for you in the file `genomics-4.nf`. Exactly as for the [nf-test side quest](../../side_quests/nf-test.md), we're going to add a few different types of tests to the three processes in this pipeline, as well as a workflow-level test.

### 0.1. Check the workflow runs

Before we start adding tests, let's make sure the workflow runs as expected.

Let's try running that now.

```bash
nextflow run genomics-4.nf -resume
```

This should look very familiar by now if you've been working through this training course from the start.

```console title="Output"
 N E X T F L O W   ~  version 24.10.0

Launching `genomics-4.nf` [gloomy_poincare] DSL2 - revision: 43203316e0

executor >  local (7)
[18/89dfa4] SAMTOOLS_INDEX (1)       | 3 of 3 ✔
[30/b2522b] GATK_HAPLOTYPECALLER (2) | 3 of 3 ✔
[a8/d2c189] GATK_JOINTGENOTYPING     | 1 of 1 ✔
```

Like previously, there will now be a `work` directory and a `results_genomics` directory inside your project directory. We'll actually make use of these results later on in our testing. But from now on we're going to be using the `nf-test` package to test the pipeline.

### 0.2. Initialize `nf-test`

As for the [nf-test side quest](../../side_quests/nf-test.md), we need to initialize the `nf-test` package.

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

### Takeaway

Now we're ready to start writing tests for our genomics pipeline.

### What's next?

Write basic tests that evaluate whether the process calls were successful and produced the correct outputs.

---

## 1. Test a process for success and matching outputs

We'll start by testing the `SAMTOOLS_INDEX` process, which creates index files for BAM files to enable efficient random access. This is a good first test case because:

1. It has a single, well-defined input (a BAM file)
2. It produces a predictable output (a BAI index file)
3. The output should be identical for identical inputs

### 1.1. Generate a test file stub

First, generate a test file stub:

```bash
nf-test generate process modules/samtools/index/main.nf
```

This creates a file in the same directory as `main.nf`, summarized in the terminal output as follows:

```console title="Output"
Load source file '/workspaces/training/nf4-science/genomics/modules/samtools/index/main.nf'
Wrote process test file '/workspaces/training/nf4-science/genomics/tests/modules/samtools/index/main.nf.test

SUCCESS: Generated 1 test files.
```

You can navigate to the directory in the file explorer and open the file, which should contain the following code:

```groovy title="tests/modules/samtools/index/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process SAMTOOLS_INDEX"
    script "modules/samtools/index/main.nf"
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

The starting assertions should be familiar from the [nf-test side quest](../../side_quests/nf-test.md):

- `assert process.success` states that we expect the process to run successfully and complete without any failures.
- `snapshot(process.out).match()` states that we expect the result of the run to be identical to the result obtained in a previous run (if applicable).
  We discuss this in more detail later.

Using this as a starting point, we need to add the right test inputs for the samtools index process, and any parameters if applicable.

### 1.2. Move the test file and update the script path

Before we get to work on filling out the test, we need to move the file to its definitive location. Part of the reason we added a directory for each module is that we can now ship tests in a `tests` directory co-located with each module's `main.nf` file. Let's create that directory and move the test file there.

```bash
mkdir -p modules/samtools/index/tests
mv tests/modules/samtools/index/main.nf.test modules/samtools/index/tests/
```

Now we can simplify the `script` section of the test file to a relative path:

_Before:_

```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
name "Test Process SAMTOOLS_INDEX"
script "modules/samtools/index/main.nf"
process "SAMTOOLS_INDEX"
```

_After:_

```groovy title="modules/samtools/index/tests/main.nf.test" linenums="3" hl_lines="2"
name "Test Process SAMTOOLS_INDEX"
script "../main.nf"
process "SAMTOOLS_INDEX"
```

This tells the test where to find the module's `main.nf` file, without having to specify the full path.

### 1.3. Provide test inputs for SAMTOOLS_INDEX

The stub file includes a placeholder that we need to replace with an actual test input, appropriate to the input of `samtools index`. The appropriate input is a BAM file, which we have available in the `data/bam` directory.

_Before:_

```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
process {
    """
    // define inputs of the process here. Example:
    // input[0] = file("test-file.txt")
    """
}
```

_After:_

```groovy title="modules/samtools/index/tests/main.nf.test" linenums="14"
process {
    """
    input[0] = file("${projectDir}/data/bam/reads_son.bam")
    """
}
```

### 1.4. Name the test based on functionality

As we learned before, it's good practice to rename the test to something that makes sense in the context of the test.

_Before:_

```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
test("Should run without failures") {
```

This takes an arbitrary string, so we could put anything we want.
Here we choose to refer to the file name and its format:

_After:_

```groovy title="modules/samtools/index/tests/main.nf.test" linenums="7"
test("Should index reads_son.bam correctly") {
```

### 1.5. Specify test parameters

The `params` block in the stub file includes a placeholder for parameters:

_Before:_

```groovy title="modules/samtools/index/tests/main.nf.test" linenums="11"
params {
    // define parameters here. Example:
    // outdir = "tests/results"
}
```

We use it to specify a location for the results to be output, using the default suggestion:

_After:_

```groovy title="modules/samtools/index/tests/main.nf.test" linenums="11"
params {
    outdir = "tests/results"
}
```

### 1.6. Run the test and examine the output

Run the test:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

This should produce:

```console title="Output"
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process SAMTOOLS_INDEX

  Test [e761f2e2] 'Should index reads_son.bam correctly' PASSED (11.226s)
  Snapshots:
    1 created [Should index reads_son.bam correctly]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 11.265s
```

As we learned previously, this verified the basic assertion about the success of the process and created a snapshot file based on the output of the process. We can see the contents of the snapshot file in the `tests/modules/samtools/index/tests/main.nf.test.snap` file:

```json title="tests/modules/samtools/index/tests/main.nf.test.snap" linenums="1"
{
  "Should index reads_son.bam correctly": {
    "content": [
      {
        "0": [
          [
            "reads_son.bam:md5,af5956d9388ba017944bef276b71d809",
            "reads_son.bam.bai:md5,a2ca7b84998218ee77eff14af8eb8ca2"
          ]
        ]
      }
    ],
    "meta": {
      "nf-test": "0.9.2",
      "nextflow": "24.10.0"
    },
    "timestamp": "2025-03-03T16:59:54.195992321"
  }
}
```

We can also run the test again and see that it passes, because the output is identical to the snapshot:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

This produces:

```console title="Output"
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process SAMTOOLS_INDEX

  Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (11.91s)


SUCCESS: Executed 1 tests in 11.947s
```

### 1.7. Add more tests to `SAMTOOLS_INDEX`

Sometimes it's useful to test a range of different input files to ensure we're testing for a variety of potential issues. Let's also test for the mother and father's bam files in the trio from our test data. Add the following tests to the test file:

```
    test("Should index reads_mother.bam correctly") {

        when {
            params {
                outdir = "tests/results"
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_mother.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }

    test("Should index reads_father.bam correctly") {

        when {
            params {
                outdir = "tests/results"
            }
            process {
                """
                input[0] = file("${projectDir}/data/bam/reads_father.bam")
                """
            }
        }

        then {
            assert process.success
            assert snapshot(process.out).match()
        }

    }
```

Then you can run the test again:

```bash
nf-test test modules/samtools/index/tests/main.nf.test
```

This produces:

```groovy title="modules/samtools/index/tests/main.nf.test" linenums="27"
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr

Warning: every snapshot that fails during this test run is re-record.

Test Process SAMTOOLS_INDEX

  Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (10.916s)
  Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (11.147s)
  Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (10.785s)
  Snapshots:
    2 created [Should index reads_father.bam correctly, Should index reads_mother.bam correctly]


Snapshot Summary:
  2 created

SUCCESS: Executed 3 tests in 32.913s
```

Notice the warning, referring to the effect of the `--update-snapshot` parameter.

!!!note

    Here we are using test data that we used previously to demonstrate the scientific outputs of the pipeline.
    If we had been planning to operate these tests in a production environment, we would have generated smaller inputs for testing purposes.

    In general it's important to keep unit tests as light as possible by using the smallest pieces of data necessary and sufficient for evaluating process functionality, otherwise the total runtime can add up quite seriously.
    A test suite that takes too long to run regularly is a test suite that's likely to get skipped in the interest of expediency.

### Takeaway

You've written your first module test for a genomics process, verifying that `SAMTOOLS_INDEX` correctly creates index files for different BAM files. The test suite ensures that:

1. The process runs successfully
2. Index files are created
3. The outputs are consistent across runs
4. The process works for all sample BAM files

### What's next?

Learn how to write tests for other processes in our genomics workflow, using the setup method to handle chained processes. We'll also evaluate whether outputs, specifically our VCF files, contain expected variant calls.

---

## 2. Add tests to a chained process and test for contents

To test `GATK_HAPLOTYPECALLER`, we need to provide the process with the `SAMTOOLS_INDEX` output as an input. We could do that by running `SAMTOOLS_INDEX`, retrieving its outputs, and storing them with the test data for the workflow. That's actually the recommended approach for a polished pipeline, but nf-test provides an alternative approach, using the `setup` method.

With the setup method, we can trigger the `SAMTOOLS_INDEX` process as part of the test setup, and then use its output as an input for `GATK_HAPLOTYPECALLER`. This has a cost - we're going to have to run the `SAMTOOLS_INDEX` process every time we run the test for `GATK_HAPLOTYPECALLER`- but maybe we're still developing the workflow and don't want to pre-generate test data we might have to change later. `SAMTOOLS_INDEX` process is also very quick, so maybe the benefits of pre-generating and storing its outputs are negligible. Let's see the setup method works.

### 2.1. Generate and place the test file

As previously, first we generate the file stub:

```bash
nf-test generate process modules/gatk/haplotypecaller/main.nf
```

This produces the following test stub:

```groovy title="tests/modules/gatk/haplotypecaller/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_HAPLOTYPECALLER"
    script "modules/gatk/haplotypecaller/main.nf"
    process "GATK_HAPLOTYPECALLER"

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

### 2.2. Move the test file and update the script path

We create a directory for the test file co-located with the module's `main.nf` file:

```bash
mkdir -p modules/gatk/haplotypecaller/tests
```

And move the test stub file there:

```bash
mv tests/modules/gatk/haplotypecaller/main.nf.test modules/gatk/haplotypecaller/tests/
```

Finally, don't forget to update the script path:

_Before:_

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_HAPLOTYPECALLER"
    script "modules/gatk/haplotypecaller/main.nf"
    process "GATK_HAPLOTYPECALLER"
```

_After:_

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="3" hl_lines="2"
    name "Test Process GATK_HAPLOTYPECALLER"
    script "../main.nf"
    process "GATK_HAPLOTYPECALLER"
```

### 2.3. Provide inputs using the setup method

We insert a `setup` block before the `when` block, where we can trigger a run of the `SAMTOOLS_INDEX` process on one of our original input files. Also, remember as before to change the test name to something meaningful.

_Before:_

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7"  hl_lines="1"
test("Should run without failures") {

    when {
```

_After:_

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="7" hl_lines="1-12"
    test("Should call son's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_son.bam")
                    """
                }
            }
        }

        when {
```

Then we can refer to the output of that process in the `when` block where we specify the test inputs:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="20"
        when {
            params {
                outdir = "tests/results"
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }
```

Make that change and run the test again:

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

This produces the following output:

```console title="Output"
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process GATK_HAPLOTYPECALLER

  Test [30247349] 'Should call son's haplotype correctly' PASSED (20.002s)
  Snapshots:
    1 created [Should call son's haplotype correctly]


Snapshot Summary:
  1 created

SUCCESS: Executed 1 tests in 20.027s
```

It also produces a snapshot file like earlier.

### 2.4. Run again and observe failure

Interestingly, if you run the exact same command again, this time the test will fail. This command:

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

produces:

```console title="Output"
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process GATK_HAPLOTYPECALLER

  Test [c847bfae] 'Should call son's haplotype correctly' FAILED (19.979s)

  java.lang.RuntimeException: Different Snapshot:
  [                                                                                           [
      {                                                                                           {
          "0": [                                                                                      "0": [
              "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
          ],                                                                                          ],
          "1": [                                                                                      "1": [
              "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
          ],                                                                                          ],
          "idx": [                                                                                    "idx": [
              "reads_son.bam.g.vcf.idx:md5,dc36c18f2afdc546f41e68b2687e9334"                     |                 "reads_son.bam.g.vcf.idx:md5,dbad4b76a4b90c158ffc9c9740764242"
          ],                                                                                          ],
          "vcf": [                                                                                    "vcf": [
              "reads_son.bam.g.vcf:md5,069316cdd4328542ffc6ae247b1dac39"                         |                 "reads_son.bam.g.vcf:md5,005f1a13ee39f11b0fc9bea094850eac"
          ]                                                                                           ]
      }                                                                                           }
  ]                                                                                           ]

  Nextflow stdout:

  Nextflow stderr:


    Obsolete snapshots can only be checked if all tests of a file are executed successful.


FAILURE: Executed 1 tests in 20.032s (1 failed)
```

The error message tells you there were differences between the snapshots for the two runs; specifically, the md5sum values are different for the VCF files.

Why? To make a long story short, the HaplotypeCaller tool includes a timestamp in the VCF header that is different every time (by definition).
As a result, we can't just expect the files to have identical md5sums even if they have identical content in terms of the variant calls themselves.

How do we deal with that?

### 2.5. Use a content assertion method to check a specific variant

One way to solve the problem is to use a [different kind of assertion](https://nf-co.re/docs/contributing/tutorials/nf-test_assertions).
In this case, we're going to check for specific content instead of asserting identity.
More exactly, we'll have the tool read the lines of the VCF file and check for the existence of specific lines.

In practice, we replace the second assertion in the `then` block as follows:

_Before:_

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3"
then {
    assert process.success
    assert snapshot(process.out).match()
}
```

_After:_

```console title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="35" hl_lines="3 4"
        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_son')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3282	GT:DP:GQ:MIN_DP:PL	0/0:25:72:24:0,72,719')
        }
```

Here we're reading in the full content of the VCF output file and searching for a content match, which is okay to do on a small test file, but you wouldn't want to do that on a larger file.
You might instead choose to read in specific lines.

This approach does require choosing more carefully what we want to use as the 'signal' to test for.
On the bright side, it can be used to test with great precision whether an analysis tool can consistently identify 'difficult' features (such as rare variants) as it undergoes further development.

### 2.6. Run again and observe success

Once we've modified the test in this way, we can run the test multiple times, and it will consistently pass.

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

This produces:

```console title="Output"
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process GATK_HAPLOTYPECALLER

  Test [c847bfae] 'Should call son's haplotype correctly' PASSED (19.33s)


SUCCESS: Executed 1 tests in 19.382s
```

### 2.7. Add more tests

Add similar tests for the mother and father samples:

```groovy title="modules/gatk/haplotypecaller/tests/main.nf.test" linenums="43"
    test("Should call mother's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_mother.bam")
                    """
                }
            }
        }

        when {
            params {
                outdir = "tests/results"
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_mother')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3278	GT:DP:GQ:MIN_DP:PL	0/0:38:99:37:0,102,1530')
        }
    }

    test("Should call father's haplotype correctly") {

        setup {
            run("SAMTOOLS_INDEX") {
                script "../../../samtools/index/main.nf"
                process {
                    """
                    input[0] =  file("${projectDir}/data/bam/reads_father.bam")
                    """
                }
            }
        }

        when {
            params {
                outdir = "tests/results"
            }
            process {
                """
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("${projectDir}/data/ref/ref.fasta")
                input[2] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[3] = file("${projectDir}/data/ref/ref.dict")
                input[4] = file("${projectDir}/data/ref/intervals.bed")
                """
            }
        }

        then {
            assert process.success
            assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father')
            assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3277	.	G	<NON_REF>	.	.	END=3281	GT:DP:GQ:MIN_DP:PL	0/0:44:99:42:0,120,1800')
        }
    }
```

### 2.8. Run the test command

```bash
nf-test test modules/gatk/haplotypecaller/tests/main.nf.test
```

This produces:

```console title="Output"
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process GATK_HAPLOTYPECALLER

  Test [c847bfae] 'Should call son's haplotype correctly' PASSED (19.91s)
  Test [44494e9c] 'Should call mother's haplotype correctly' PASSED (18.606s)
  Test [eb0d1a07] 'Should call father's haplotype correctly' PASSED (18.773s)


SUCCESS: Executed 3 tests in 57.348s
```

That completes the basic test plan for this second step in the pipeline. On to the third and last module-level test!

### Takeaway

You've learned how to:

1. Test processes that depend on outputs from other processes
2. Verify specific genomic variants in VCF output files
3. Handle non-deterministic outputs by checking specific content
4. Test variant calling across multiple samples

### What's next?

Learn how to write tests that use pre-generated test data for the joint genotyping step.

---

## 3. Use pre-generated test data

For the joint genotyping step, we'll use a different approach - using pre-generated test data. This is often preferable for:

1. Complex processes with multiple dependencies
2. Processes that take a long time to run
3. Processes that are part of a stable, production pipeline

### 3.1. Generate test data

Let's inspect the results we generated at the start of this section:

```bash
tree results_genomics/
```

```console title="Results directory contents"
results_genomics/
├── family_trio.joint.vcf
├── family_trio.joint.vcf.idx
├── reads_father.bam -> /workspaces/training/nf4-science/genomics/work/42/a3bf19dbfaf1f3672b16a5d5e6a8be/reads_father.bam
├── reads_father.bam.bai -> /workspaces/training/nf4-science/genomics/work/cf/289c2d264f496d60a69e3e9ba6463e/reads_father.bam.bai
├── reads_father.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf
├── reads_father.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/30/b2522b83c63baff8c3cf75704512a2/reads_father.bam.g.vcf.idx
├── reads_mother.bam -> /workspaces/training/nf4-science/genomics/work/af/f31a6ade82cc0cf853c4f61c8bc473/reads_mother.bam
├── reads_mother.bam.bai -> /workspaces/training/nf4-science/genomics/work/18/89dfa40a3def17e45421e54431a126/reads_mother.bam.bai
├── reads_mother.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf
├── reads_mother.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/f6/be2efa58e625d08cf8d0da1d0e9f09/reads_mother.bam.g.vcf.idx
├── reads_son.bam -> /workspaces/training/nf4-science/genomics/work/9f/9615dd553d6f13d8bec4f006ac395f/reads_son.bam
├── reads_son.bam.bai -> /workspaces/training/nf4-science/genomics/work/4d/cb384a97db5687cc9daab002017c7c/reads_son.bam.bai
├── reads_son.bam.g.vcf -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf
└── reads_son.bam.g.vcf.idx -> /workspaces/training/nf4-science/genomics/work/fe/2f22d56aa16ed45f8bc419312894f6/reads_son.bam.g.vcf.idx

0 directories, 14 files
```

Note that some of these files are symlinks to the actual files in the `work` directory.

The joint genotyping step needs the VCF files produced by the haplotype caller steps as inputs, along with the indices. So let's copy the results we have into the `jointgenotyping` module's test directory. We'll follow symlinks to the original files.

```bash
mkdir -p modules/gatk/jointgenotyping/tests/inputs/
cp -rL results_genomics/*.g.vcf results_genomics/*.g.vcf.idx modules/gatk/jointgenotyping/tests/inputs/
```

Now we can use these files as inputs to the test we're going to write for the joint genotyping step.

### 3.2. Generate the test file stub

As previously, first we generate the file stub:

```bash
nf-test generate process modules/gatk/jointgenotyping/main.nf
```

This produces the following test stub:

```groovy title="tests/modules/gatk/jointgenotyping/main.nf.test" linenums="1"
nextflow_process {

    name "Test Process GATK_JOINTGENOTYPING"
    script "modules/gatk/jointgenotyping/main.nf"
    process "GATK_JOINTGENOTYPING"

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

### 3.3. Move the test file and update the script path

This time we already have a directory for tests co-located with the module's `main.nf` file, so we can move the test stub file there:

```bash
mv tests/modules/gatk/jointgenotyping/main.nf.test modules/gatk/jointgenotyping/tests/
```

And don't forget to update the script path:

_Before:_

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
name "Test Process GATK_JOINTGENOTYPING"
script "modules/gatk/jointgenotyping/main.nf"
process "GATK_JOINTGENOTYPING"
```

_After:_

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="3" hl_lines="2"
name "Test Process GATK_JOINTGENOTYPING"
script "../main.nf"
process "GATK_JOINTGENOTYPING"
```

### 3.4. Provide inputs

Fill in the inputs based on the process input definitions and rename the test accordingly:

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="7"
    test("Should call trio's joint genotype correctly") {

        when {
            params {
                outdir = "tests/results"
            }
            process {
                """
                input[0] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf")
                ]
                input[1] = [
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_father.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_mother.bam.g.vcf.idx"),
                    file("${projectDir}/modules/gatk/jointgenotyping/tests/inputs/reads_son.bam.g.vcf.idx")
                ]
                input[2] = file("${projectDir}/data/ref/intervals.bed")
                input[3] = "family_trio"
                input[4] = file("${projectDir}/data/ref/ref.fasta")
                input[5] = file("${projectDir}/data/ref/ref.fasta.fai")
                input[6] = file("${projectDir}/data/ref/ref.dict")
                """
            }
        }
```

### 3.5. Use content assertions

The output of the joint genotyping step is another VCF file, so we're going to use a content assertion again.

```groovy title="modules/gatk/jointgenotyping/tests/main.nf.test" linenums="25"
    then {
        assert process.success
        assert path(process.out[0][0]).readLines().contains('#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	reads_father	reads_mother	reads_son')
        assert path(process.out[0][0]).readLines().contains('20_10037292_10066351	3480	.	C	CT	1625.89	.	AC=5;AF=0.833;AN=6;BaseQRankSum=0.220;DP=85;ExcessHet=0.0000;FS=2.476;MLEAC=5;MLEAF=0.833;MQ=60.00;MQRankSum=0.00;QD=21.68;ReadPosRankSum=-1.147e+00;SOR=0.487	GT:AD:DP:GQ:PL	0/1:15,16:31:99:367,0,375	1/1:0,18:18:54:517,54,0	1/1:0,26:26:78:756,78,0')
    }
```

By checking the content of a specific variant in the output file, this test verifies that:

1. The joint genotyping process runs successfully
2. The output VCF contains all three samples in the correct order
3. A specific variant is called correctly with:
   - Accurate genotypes for each sample (0/1 for father, 1/1 for mother and son)
   - Correct read depths and genotype qualities
   - Population-level statistics like allele frequency (AF=0.833)

We haven't snapshotted the whole file, but by checking a specific variant, we can be confident that the joint genotyping process is working as expected.

### 3.6. Run the test

```bash
nf-test test modules/gatk/jointgenotyping/tests/main.nf.test
```

This produces:

```console title="Output"
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process GATK_JOINTGENOTYPING

  Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (21.604s)


SUCCESS: Executed 1 tests in 21.622s
```

The test passes, verifying that our joint genotyping process correctly:

1. Combines individual sample VCFs
2. Performs joint variant calling
3. Produces a multi-sample VCF with consistent genotype calls across runs

### Takeaway

You know how to:

- Use previously generated results as inputs for tests
- Write tests using pre-generated test data

### What's next?

Add a workflow-level test to verify the entire variant calling pipeline works end-to-end.

---

## 4. Add a workflow-level test

Now we'll test the complete variant calling pipeline, from BAM files to joint genotypes. This verifies that:

1. All processes work together correctly
2. Data flows properly between steps
3. Final variant calls are consistent

### 4.1. Generate the workflow test

Generate a test file for the complete pipeline:

```bash
nf-test generate pipeline genomics-4.nf
```

This creates a basic test stub:

```groovy title="tests/genomics-4.nf.test" linenums="1"
nextflow_pipeline {

    name "Test Workflow genomics-4.nf"
    script "genomics-4.nf"

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

Just correct the name to something meaningful (you'll see why this is useful shortly).

_Before:_

```groovy title="tests/genomics-4.nf.test" linenums="1" hl_lines="1"
    test("Should run without failures") {
```

_After:_

```groovy title="tests/genomics-4.nf.test" linenums="1" hl_lines="1"
    test("Should run the pipeline without failures") {
```

!!!note

    In this case the test file can stay where `nf-test` created it.

### 4.2. Specify input parameters

We still need to specify inputs, which is done slightly different at the workflow level compared to module-level tests.
There are several ways of doing this, including by specifying a profile.
However, a simpler way is to set up a `params {}` block in the `nextflow.config` file that `nf-test init` originally created in the `tests` directory.

```groovy title="tests/nextflow.config" linenums="1"
/*
 * Pipeline parameters
 */

params {
    // Primary input (file of input files, one per line)
    reads_bam        = "${projectDir}/data/sample_bams.txt"

    // Output directory
    outdir = "results_genomics"

    // Reference genome and intervals
    reference        = "${projectDir}/data/ref/ref.fasta"
    reference_index  = "${projectDir}/data/ref/ref.fasta.fai"
    reference_dict   = "${projectDir}/data/ref/ref.dict"
    intervals        = "${projectDir}/data/ref/intervals.bed"

    // Base name for final output file
    cohort_name      = "family_trio"
}
```

When we run the test, `nf-test` will pick up this configuration file and pull in the inputs accordingly.

### 4.3. Run the workflow test

```bash
nf-test test tests/genomics-4.nf.test
```

This produces:

```console title="Output"
🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Workflow genomics-4.nf

  Test [c7dbcaca] 'Should run without failures' PASSED (48.443s)


SUCCESS: Executed 1 tests in 48.486s
```

The test passes, confirming that our complete variant calling pipeline:

1. Successfully processes all samples
2. Correctly chains together all steps

### 4.4. Run ALL tests

nf-test has one more trick up it's sleeve. We can run all the tests at once! Modify the `nf-test.config` file so that nf-test looks in every directory for nf-test files. You can do this by modifying the `testsDir` parameter:

_Before:_

```groovy title="nf-test.config" linenums="1" hl_lines="3"
config {

    testsDir "tests"
    workDir ".nf-test"
    configFile "tests/nextflow.config"
    profile ""

}
```

_After:_

```groovy title="nf-test.config" linenums="1" hl_lines="3"
config {

    testsDir "."
    workDir ".nf-test"
    configFile "tests/nextflow.config"
    profile ""

}
```

Now, we can simply run nf-test and it will run _every single test_ in our repository:

```bash
nf-test test
```

This produces:

```console title="Output"

🚀 nf-test 0.9.2
https://www.nf-test.com
(c) 2021 - 2024 Lukas Forer and Sebastian Schoenherr


Test Process GATK_HAPLOTYPECALLER

  Test [c847bfae] 'Should call son's haplotype correctly' PASSED (20.951s)
  Test [44494e9c] 'Should call mother's haplotype correctly' PASSED (19.155s)
  Test [eb0d1a07] 'Should call father's haplotype correctly' PASSED (21.843s)

Test Process GATK_JOINTGENOTYPING

  Test [ac2067de] 'Should call trio's joint genotype correctly' PASSED (22.994s)

Test Process SAMTOOLS_INDEX

  Test [625e39ee] 'Should index reads_son.bam correctly' PASSED (11.281s)
  Test [a8b28f36] 'Should index reads_mother.bam correctly' PASSED (11.126s)
  Test [c15852a1] 'Should index reads_father.bam correctly' PASSED (12.005s)

Test Workflow genomics-4.nf

  Test [c7dbcaca] 'Should run the pipeline without failures' PASSED (47.92s)


SUCCESS: Executed 8 tests in 167.772s
```

7 tests in 1 command! We spent a long time configuring lots and lots of tests, but when it came to running them it was very quick and easy. You can see how useful this is when maintaining a large pipeline, which could include hundreds of different elements. We spend time writing tests once so we can save time running them many times.

Furthermore, we can automate this! Imagine tests running every time you or a colleague tries to add new code. This is how we ensure our pipelines maintain a high standard.

## Takeaway

You now know how to write and run several kinds of tests for your genomics pipeline using nf-test. This testing framework helps ensure your variant calling workflow produces consistent, reliable results across different environments and as you make code changes.

You've learned to test critical components like:

- The `SAMTOOLS_INDEX` process that prepares BAM files for variant calling
- The `GATK_HAPLOTYPECALLER` process that identifies variants in individual samples
- The `GATK_JOINTGENOTYPING` process that combines variant calls across a cohort

You've also implemented different testing strategies specific to genomics data:

- Verifying that VCF files contain expected variant calls despite non-deterministic elements like timestamps
- Testing with a family trio dataset to ensure proper variant identification across related samples
- Checking for specific genomic coordinates and variant information in your output files

These testing skills are essential for developing robust bioinformatics pipelines that can reliably process genomic data and produce accurate variant calls. As you continue working with Nextflow for genomics analysis, this testing foundation will help you maintain high-quality code that produces trustworthy scientific results.
