# Part 3: Run pipelines from remote repositories

So far, you've run workflow scripts stored locally.
In practice, you'll often want to run pipelines published in remote repositories, such as GitHub, without downloading them yourself.

Nextflow makes this straightforward: you can run any pipeline directly from a Git repository URL.

---

## 1. Run a pipeline from GitHub

The basic syntax for running a remote pipeline is `nextflow run <repository>`, where `<repository>` can be a GitHub repository path like `nextflow-io/hello`, a full URL, or a path to GitLab, Bitbucket, or another Git hosting service.

### 1.1. Launch the pipeline

Run the official Nextflow "hello" demo pipeline.

```bash
nextflow run nextflow-io/hello
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Pulling nextflow-io/hello ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sleepy_swanson] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

The first time you run a remote pipeline, Nextflow downloads it and caches it locally.
Subsequent runs reuse the cached version unless you explicitly request an update.

### Takeaway

You know how to run a pipeline directly from a GitHub repository without downloading it yourself.

### What's next?

Learn how to pin a specific version of a remote pipeline for reproducibility.

---

## 2. Specify a version for reproducibility

By default, Nextflow runs the latest revision from the default branch.
You can pin a particular version (tag), branch, or commit using the `-r` flag.

### 2.1. Pin a specific revision

```bash
nextflow run nextflow-io/hello -r v1.3
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Pinning an exact revision is essential for reproducibility.
It guarantees that you and your collaborators run the exact same pipeline code, regardless of what has changed in the repository since.

### Takeaway

You know how to pin a remote pipeline to a specific version, branch, or commit for reproducible execution.

### What's next?

You've covered the fundamentals of running Nextflow pipelines.
See [Next steps](next_steps.md) for where to go from here.

---

## Summary

In this part you learned to:

- Run a pipeline directly from a GitHub repository without downloading it
- Pin a remote pipeline to a specific revision for reproducibility
