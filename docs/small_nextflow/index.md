---
title: Small Nextflow
hide:
  - toc
---

# Small Nextflow

Hello! Welcome to Small Nextflow, a hands-on workshop where you'll build a real-world image classification workflow from scratch.

In this workshop, you'll create a workflow that fetches cat images, classifies them using machine learning, and produces visual collages.
Along the way, you'll learn the core concepts that make Nextflow powerful for scientific computing: channels, processes, operators, and reproducible execution.

By starting with an empty directory and building up piece by piece, you'll gain an intuitive understanding of how Nextflow workflows are structured and how data flows through them.

Let's get started!

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

## Learning objectives

In this workshop, you will learn foundational Nextflow concepts by building a complete workflow from scratch.

By the end of this workshop you will be able to:

- Create channels from files and understand channel types (queue vs. value channels)
- Define processes with input, output, and script blocks
- Use operators to transform and combine channel data (`map`, `join`, `groupTuple`, `collect`)
- Work with metadata using tuples and maps
- Configure workflows with parameters
- Manage resources and direct process execution
- Publish workflow outputs in organized structures
- Containerize workflows for reproducibility
- Make workflows portable across filesystems

## Audience & prerequisites

This workshop is designed for those who want to learn Nextflow by building a complete workflow from the ground up.
Some basic familiarity with the command line and programming concepts is helpful but not required.

**Prerequisites**

- A GitHub account
- Basic familiarity with command line
- Curiosity about machine learning (no ML expertise needed!)

## Workshop structure

The workshop is organized into four chapters:

**Chapter 1: Fundamentals** - Learn the building blocks of Nextflow by creating channels, defining processes, and working with parameters and metadata.

**Chapter 2: Data Transformation & Analysis** - Build a multi-step workflow that classifies images, manages resources, and groups results intelligently.

**Chapter 3: Publishing & Portability** - Make your workflow production-ready by publishing organized outputs, supporting multiple filesystems, and containerizing dependencies.

**Chapter 4: Advanced Topics** - Explore version control integration, cloud execution, and extension exercises to deepen your skills.
