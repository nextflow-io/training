# Orientation

The Gitpod environment contains some test data that will be used in this training course. All software required is already installed and configured in it too.

!!! note

    Follow [this link](../../envsetup/) if you have not yet set up your Gitpod environment.

## Materials provided

Throughout this training course, we'll be working in the `hello-nextflow/` directory.

```bash
cd /workspace/gitpod/hello-nextflow
```

This directory contains all the code files, test data and accessory files you will need. Feel free to explore the contents of this directory; an easy way to see what it contains is the use the `tree` command (here we generate a table of contents to the second level down).

```bash
tree . -L 2
```

You should see the following output:

```console title="Directory contents"
/workspace/gitpod/hello-nextflow
├── data
│   ├── bam
│   ├── greetings.csv
│   ├── ref
│   ├── sample_bams.txt
│   └── samplesheet.csv
├── hello-gatk.nf
├── hello-modules.nf
├── hello-nf-test.nf
├── hello-world.nf
├── nextflow.config
└── scripts
    ├── hello-config-1.config
    ├── hello-config-2.config
    ├── hello-config-3.config
    ├── hello-config-4.config
    ├── hello-gatk-1.nf
    ├── hello-gatk-2.nf
    ├── hello-gatk-3.nf
    ├── hello-gatk-4.nf
    ├── hello-gatk-5.nf
    ├── hello-gatk-6.nf
    ├── hello-modules-1.nf
    ├── hello-modules-2.nf
    ├── hello-modules-3.nf
    ├── hello-world-1.nf
    ├── hello-world-2.nf
    ├── hello-world-3.nf
    ├── hello-world-4.nf
    ├── hello-world-5.nf
    ├── hello-world-6.nf
    ├── hello-world-7.nf
    ├── hello-world-8.nf
    ├── hello-world-9.nf
    ├── modules
    └── nextflow.config

13 directories, 48 files

```

**The `data` directory** contains the input data we'll use in Part 2: Hello Science, which uses an example from genomics to demonstrate how to build a simple analysis pipeline. The data is described in detail in that section of the course.

**The file `nextflow.config`** is a configuration file that sets minimal environment properties.

**The file `hello-world.nf`** is a simple but fully functional workflow script that serves as a starting point to Part 1: Hello World.

**The file `hello-gatk.nf`** is a stub that serves as a starting point to Part 2: Hello Science. In its initial state, it is NOT a functional workflow script.

**The remaining `.nf` files** are functional workflow scripts that serve as starting points for the corresponding parts of the course.

**The `scripts` directory** contains the completed workflow scripts that result from each step of the course. They are intended to be used as a reference to check your work and troubleshoot any issues. The name and number in the filename correspond to the step of the relevant part of the course. For example, the file `hello-world-4.nf` is the expected result of completing steps 1 through 4 of Part 1: Hello World.
