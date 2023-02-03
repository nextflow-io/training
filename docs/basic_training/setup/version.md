---
title: Nextflow version
description: How to use a specific version of Nextflow
---

# Selecting a Nextflow version

By default, Nextflow will pull the latest stable version into your environment.
Nextflow is constantly evolving as we make improvements and fix bugs.
The latest releases can be viewed on GitHub: <https://github.com/nextflow-io/nextflow>

If you want to use a specific version of Nextflow, you can set the `NXF_VER` variable as shown below:

```bash
export NXF_VER=22.04.5
```

!!! Note

    Most of this tutorial workshop requires `NXF_VER=22.04.0` or later, to use DSL2 as default.

Run `nextflow -version` again to confirm that the change has taken effect.
