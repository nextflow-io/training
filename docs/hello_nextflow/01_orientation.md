# Orientation

The Gitpod environment contains some test data that will be used in this workshop. All software required are already installed and configured in it too.

!!! note

    Follow [this link](../../envsetup/) if you have not yet setup your Gitpod environment.

## Getting started

You will complete this module in the `hello-nextflow/` folder.

```bash
cd /workspace/gitpod/hello-nextflow
```

In this folder you will all test data, code and accessory needed to work through this training module.

!!! question "Exercise"

    View all the folder and files in the `hello-nextflow` directory.

    ```console
    tree .
    ```

You should see the following output:

```console title="Output"
/workspace/gitpod/hello-nextflow
├── data
│   ├── bam
│   │   ├── reads_father.bam
│   │   ├── reads_mother.bam
│   │   └── reads_son.bam
│   ├── greetings.txt
│   ├── intervals.list
│   ├── ref.tar.gz
│   ├── sample_bams.txt
│   └── samplesheet.csv
├── hello-gatk.nf
├── hello-modules.nf
├── hello-nf-test.nf
├── hello-world.nf
├── nextflow.config
└── scripts
    ├── hello-gatk-1.nf
    ├── hello-gatk-2.nf
    ├── hello-gatk-3.nf
    ├── hello-gatk-4.nf
    ├── hello-gatk-5.nf
    ├── hello-gatk-6.nf
    ├── hello-modules-1.nf
    ├── hello-modules-2.nf
    ├── hello-modules-3.nf
    ├── hello-world-10.nf
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
    │   └── local
    │       ├── gatk
    │       │   ├── haplotypecaller
    │       │   │   └── main.nf
    │       │   └── jointgenotyping
    │       │       ├── main.nf
    │       │       └── tests
    │       │           └── inputs
    │       │               ├── family_trio_map.tsv
    │       │               ├── reads_father.bam.g.vcf
    │       │               ├── reads_father.bam.g.vcf.idx
    │       │               ├── reads_mother.bam.g.vcf
    │       │               ├── reads_mother.bam.g.vcf.idx
    │       │               ├── reads_son.bam.g.vcf
    │       │               └── reads_son.bam.g.vcf.idx
    │       └── samtools
    │           └── index
    │               └── main.nf
    └── nextflow.config

12 directories, 43 files

```

Each file will be used in this training module.

**The `data` directory** contains the input data we'll use in Part 2: Hello GATK, which uses an example from genomics to demonstrate how to build a simple analysis pipeline. The data is described in detail in that section of the training.

**The `scripts` directory** contains the completed workflow scripts that result from each step of the tutorial and are intended to be used as a reference to check your work. The name and number in the filename correspond to the step of the relevant tutorial. For example, the file `hello-world-4.nf` is the expected result of completing steps 1 through 4 of Part 1: Hello World.

**The file `greetings.txt`** is a plain text file used to provide inputs in Part 1: Hello World.

**The file `hello-gatk.nf`** is a stub that serves as a starting point to Part 2: Hello GATK. In its initial state, it is NOT a functional workflow script.

**The file `hello-world.nf`** is a simple but fully functional workflow script that serves as a starting point to Part 1: Hello World.

**The file `nextflow.config`** is a configuration file that sets minimal environment properties.
