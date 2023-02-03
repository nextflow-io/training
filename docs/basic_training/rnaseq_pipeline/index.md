---
title: Introduction
description: Basic Nextflow Training Workshop
---

# Simple RNA-Seq pipeline

To demonstrate a real-world biomedical scenario, we will implement a proof of concept RNA-Seq pipeline which:

1. Indexes a transcriptome file
2. Performs quality controls
3. Performs quantification
4. Creates a MultiQC report

This will be done using a series of seven scripts, each of which builds on the previous to create a complete workflow. You can find these in the tutorial folder (`script1.nf` - `script7.nf`).
