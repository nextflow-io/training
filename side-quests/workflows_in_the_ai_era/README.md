# Workflows in the AI era

This side quest teaches why workflow managers still matter when an AI agent can run your analysis or generate the pipeline for you.

## Structure

```
workflows_in_the_ai_era/
├── bash/                    # Part 1: build scripts here
├── data/
│   └── samples.csv          # sample metadata (real RNA-seq data URLs)
└── nextflow/                # Part 2: build workflow here
    ├── main.nf              # starter template
    ├── modules/             # add process modules here
    └── nextflow.config      # starter configuration
```

## Tutorial parts

### Part 1: Building a bash pipeline

You'll build an RNA-seq pipeline from scratch in bash, experiencing:

- Tool installation friction (conda/mamba setup)
- Sequential processing limitations
- Parallelisation complexity (`&` and `wait`)
- Failure recovery challenges
- Reproducibility problems
- Portability issues

### Part 2: Building a Nextflow pipeline

You'll rebuild the same pipeline in Nextflow, seeing how each problem is solved:

- Container-based tool management
- Automatic parallelisation from data flow
- Built-in resume capability
- Declarative resource management
- Configuration profiles for any platform
- The workflow output system for stable result layouts

## Solutions

Complete solutions are in `../solutions/workflows_in_the_ai_era/`:

- `bash/`: final bash scripts (`process_sample.sh`, `pipeline_sequential.sh`, `pipeline_parallel.sh`)
- `nextflow/`: final Nextflow workflow with all modules

## Prerequisites

- Docker (for Nextflow containers)
- Conda or mamba (for Part 1 bash tools)
- Basic command-line familiarity

## Data

This tutorial uses real RNA-seq data from the nf-core test-datasets repository. The FASTQ files are small (about 30 seconds per sample) to keep execution fast while using actual bioinformatics tools.

## Getting started

See the full tutorial: [Workflows in the AI era](https://training.nextflow.io/side_quests/workflows_in_the_ai_era/).
