# Part 4: Run pipelines from remote repositories

So far, you've run workflow scripts stored locally.
In practice, you'll often want to run pipelines published in remote repositories, such as GitHub, without downloading them yourself.

Nextflow makes this straightforward: you can run any pipeline directly from a Git repository URL.

---

## 1. Run a pipeline from GitHub

The basic syntax for running a remote pipeline is `nextflow run <repository>`, where `<repository>` can be a GitHub repository path like `nextflow-io/hello`, a full URL, or a path to GitLab, Bitbucket, or another Git hosting service.

### 1.1. Launch the pipeline

Run the official Nextflow "hello" demo pipeline.
This is a different, much simpler pipeline than the one you've been running in this course: it predates the "Hello" pipeline used throughout this training, and just prints a greeting for each of a few hardcoded languages, so don't expect the CSV input or ASCII art you're used to.

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

### 1.2. Find where the pipeline is cached

The first time you run a remote pipeline, Nextflow downloads it and caches it locally.
Subsequent runs reuse the cached version unless you explicitly request an update.

By default, Nextflow saves pulled pipelines under `$NXF_HOME/assets`.
To find where a specific pipeline landed, and which revisions are available, ask Nextflow directly:

```bash
nextflow info nextflow-io/hello
```

??? success "Command output"

    ```console
    project name: nextflow-io/hello
    repository  : https://github.com/nextflow-io/hello
    local path  : /workspaces/.nextflow/assets/nextflow-io/hello
    main script : main.nf
    revisions   :
    * master (default)
      mybranch
      testing
      v1.1 [t]
      v1.2 [t]
      v1.3 [t]
    ```

You can also list every pipeline you've pulled so far with `nextflow list`:

```bash
nextflow list
```

??? success "Command output"

    ```console
    nextflow-io/hello
    ```

The [Use nf-core](../nfcore_run/01_run_demo.md#12-retrieve-the-pipeline-code) course covers this caching mechanism in more depth, including how to browse a pulled pipeline's source code.

### Takeaway

You know how to run a pipeline directly from a GitHub repository without downloading it yourself, and where to find it locally afterwards.

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

### 2.2. Unpin a pipeline

Pinning a revision doesn't just apply to that one run: Nextflow checks out that revision in the local cache, so it also becomes what any later run without `-r` uses.
Try running the pipeline again without `-r`:

```bash
nextflow run nextflow-io/hello
```

??? failure "Command output"

    ```console
    Project `nextflow-io/hello` is currently stuck on revision: v1.3 -- you need to explicitly specify a revision with the option `-r` in order to use it
    ```

Nextflow refuses to guess, since silently running a different revision than the one you pinned would defeat the purpose of pinning it in the first place.
To go back to running the default branch, pass it explicitly with `-r`:

```bash
nextflow run nextflow-io/hello -r master
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] DSL2 - revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

You can find the name of a pipeline's default branch by running `nextflow info <pipeline>`; it's the one marked `(default)`.
Once you've run it, the pipeline is unstuck: plain `nextflow run nextflow-io/hello` commands go back to using that branch until you pin it to something else again.

### Takeaway

You know how to pin a remote pipeline to a specific version, branch, or commit for reproducible execution, and how to unpin it again.

### What's next?

You've covered the fundamentals of running and managing Nextflow pipelines.
See [Course summary](next_steps.md) for where to go from here.

---

## Summary

In this part you learned to:

- Run a pipeline directly from a GitHub repository without downloading it
- Pin a remote pipeline to a specific revision for reproducibility, and unpin it again
