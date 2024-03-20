# Orientation

## Tour of Gitpod

If you haven't yet, log into the [![Nextflow Training GitPod](https://img.shields.io/badge/Gitpod-%20Open%20in%20Gitpod-908a85?logo=gitpod)](https://gitpod.io/#https://github.com/nextflow-io/training), which provides a virtual machine with everything already set up for you.

In the Gitpod window, you'll see a terminal. Type the following command to switch to the folder of this training material:

```bash
cd /workspace/gitpod/hello-nextflow
```

Take a few minutes to familiarize yourself with the gitpod environment, especially the file explorer, file browser and terminal.

## Pipeline data and scripts

We provide all test data, code and accessory needed to work through this training module. To view a full list, run the following command in the Gitpod terminal:

```bash
tree /workspace/gitpod/hello-nextflow
```

You should see the following output:

```console title="Output"
hello-nextflow
├── data
│   ├── bam
│       ├── reads_father.bam
│       ├── reads_mother.bam
│       └── reads_son.bam
│   ├── intervals.list
│   ├── ref.tar.gz
│   ├── sample_bams.txt
│   └── samplesheet.csv
├── scripts
│   ├── hello-gatk-1.nf
│   ├── hello-gatk-2.nf
│   ├── hello-gatk-3.nf
│   ├── hello-gatk-4.nf
│   ├── hello-gatk-5.nf
│   ├── hello-gatk-6.nf
│   ├── hello-world-1.nf
│   ├── hello-world-2.nf
│   ├── hello-world-3.nf
│   ├── hello-world-4.nf
│   ├── hello-world-5.nf
│   ├── hello-world-6.nf
│   ├── hello-world-7.nf
│   └── hello-world-8.nf
├── greetings.txt
├── hello-gatk.nf
├── hello-world.nf
└── nextflow.config

```

### Description of contents

**The `data` directory** contains the input data we'll use in Part 2: Hello GATK, which uses an example from genomics to demonstrate how to build a simple analysis pipeline. The data is described in detail in that section of the training.

**The `scripts` directory** contains the completed workflow scripts that result from each step of the tutorial and are intended to be used as a reference to check your work. The name and number in the filename correspond to the step of the relevant tutorial. For example, the file `hello-world-4.nf` is the expected result of completing steps 1 through 4 of Part 1: Hello World.

**The file `greetings.txt`** is a plain text file used to provide inputs in Part 1: Hello World.

**The file `hello-gatk.nf`** is a stub that serves as a starting point to Part 2: Hello GATK. In its initial state, it is NOT a functional workflow script.

**The file `hello-world.nf`** is a simple but fully functional workflow script that serves as a starting point to Part 1: Hello World.

**The file `nextflow.config`** is a configuration file that sets minimal environment properties.
