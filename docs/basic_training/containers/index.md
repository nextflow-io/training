---
title: Introduction
description: Basic Nextflow Training Workshop
---

# Manage dependencies and containers

Computational workflows are rarely composed of a single script or tool. Most of the time they require the usage of dozens of different software components or libraries.

Installing and maintaining such dependencies is a challenging task and the most common source of irreproducibility in scientific applications.

To overcome these issues we use containers that allow the encapsulation of software dependencies, i.e. tools and libraries required by a data analysis application in one or more self-contained, ready-to-run, immutable Linux container images, that can be easily deployed in any platform supporting the container runtime.

Containers can be executed in an isolated manner from the hosting system. Having its own copy of the file system, processing space, memory management, etc.

!!! info

    Containers were first introduced with kernel 2.6 as a Linux feature known as _Control Groups_ or [Cgroups](https://en.wikipedia.org/wiki/Cgroups).
