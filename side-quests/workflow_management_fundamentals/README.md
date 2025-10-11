# Workflow Management Fundamentals

This side quest demonstrates the benefits of workflow management systems like Nextflow by comparing a traditional bash script approach with a Nextflow workflow.

## Overview

This tutorial uses a bacterial genome analysis pipeline as an exemplar to show:

- Automatic parallelization
- Software environment isolation via containers
- Resume capability
- Portable execution across environments
- Resource management
- Data provenance tracking

## Files

- `process_samples.sh` - Traditional bash script approach
- `main.nf` - Nextflow workflow
- `modules/` - Modular process definitions
- `data/samples.csv` - Sample metadata
- `nextflow.config` - Configuration file

## Note

This is a pedagogical example. The actual sequence data files are not included as they would be too large. The focus is on understanding the workflow patterns and benefits, not executing the actual analysis.

## Getting Started

See the full tutorial in the documentation: [docs/side_quests/workflow_management_fundamentals.md](../../docs/side_quests/workflow_management_fundamentals.md)
