# Workflow Management Fundamentals

This side quest teaches why workflow managers exist by having you experience the pain points firsthand.

## Structure

```
workflow_management_fundamentals/
├── bash/                    # Part 1: Build scripts here
├── data/
│   └── samples.csv          # Sample metadata (real RNA-seq data URLs)
└── nextflow/                # Part 2: Build workflow here
    ├── main.nf              # Starter template
    ├── modules/             # Add process modules here
    └── nextflow.config      # Starter configuration
```

## Tutorial Parts

### Part 1: Building a Bash Pipeline
You'll build an RNA-seq pipeline from scratch in bash, experiencing:
- Tool installation friction (conda/mamba setup)
- Sequential processing limitations
- Parallelization complexity (& and wait)
- Failure recovery challenges
- Reproducibility problems
- Portability issues

### Part 2: Building a Nextflow Pipeline
You'll rebuild the same pipeline in Nextflow, seeing how each problem is solved:
- Container-based tool management
- Automatic parallelization from data flow
- Built-in resume capability
- Declarative resource management
- Configuration profiles for any platform
- Automatic provenance tracking

## Solutions

Complete solutions are in `../solutions/workflow_management_fundamentals/`:
- `bash/` - Final bash scripts (process_sample.sh, pipeline_sequential.sh, pipeline_parallel.sh)
- `nextflow/` - Final Nextflow workflow with all modules

## Prerequisites

- Docker (for Nextflow containers)
- Conda/Mamba (for Part 1 bash tools)
- Basic command line familiarity

## Data

This tutorial uses real RNA-seq data from the nf-core test-datasets repository.
The FASTQ files are small (~30 seconds per sample) to keep execution fast while using actual bioinformatics tools.

## Getting Started

See the full tutorial: [Workflow Management Fundamentals](https://training.nextflow.io/side_quests/workflow_management_fundamentals/)
