# Part 4: Run remote pipelines

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
     local path  : /workspaces/.nextflow/assets/.repos/nextflow-io/hello
     main script : main.nf
     revisions   :
     > master (default)
       mybranch
       testing
       v1.1 [t]
       v1.2 [t]
       v1.3 [t]
    ```

    Nextflow marks each revision you've already checked out locally with `>`; the rest are available but not yet fetched into a working copy.

You can also list every pipeline you've pulled so far with `nextflow list`:

```bash
nextflow list
```

??? success "Command output"

    ```console
    nextflow-io/hello
    ```

The [Use nf-core](../nfcore_use/01_run_demo.md#12-retrieve-the-pipeline-code) course covers this caching mechanism in more depth, including how to browse a pulled pipeline's source code.

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

    Pulling nextflow-io/hello:v1.3 ...
     downloaded from https://github.com/nextflow-io/hello.git
    Launching `https://github.com/nextflow-io/hello` [sick_carson] revision: 2ce0b0e294 [v1.3]

    executor >  local (4)
    [61/e11f77] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Bonjour world!

    Hello world!

    Hola world!
    ```

Nextflow fetches this revision the first time you request it, hence the `Pulling` and `downloaded from` lines; requesting the same revision again later skips straight to `Launching`.
Pinning an exact revision is essential for reproducibility.
It guarantees that you and your collaborators run the exact same pipeline code, regardless of what has changed in the repository since.

### 2.2. Revisions apply per invocation only

Pinning a revision with `-r` only affects the run where you specify it: it does not change what a later, plain `nextflow run` uses.
Try running the pipeline again without `-r`:

```bash
nextflow run nextflow-io/hello
```

??? success "Command output"

    ```console
    N E X T F L O W   ~  version 26.04.4

    Launching `https://github.com/nextflow-io/hello` [hungry_maxwell] revision: 3c2cdc9823 [master]

    executor >  local (4)
    [ba/08236d] sayHello (4) | 4 of 4 ✔
    Ciao world!

    Hello world!

    Bonjour world!

    Hola world!
    ```

Even though the previous run explicitly pinned `v1.3`, this run goes straight back to the default branch (`master`).
Nextflow keeps a separate local working copy for each revision you've used, which is what the `>` markers in `nextflow info` show, but it never remembers which one you ran last.
You can find the name of a pipeline's default branch by running `nextflow info <pipeline>`; it's the one marked `(default)`.
Reproducibility is entirely on you: always pass `-r` explicitly whenever it matters, rather than assuming a revision you pinned in an earlier run still applies.

### Takeaway

You know how to pin a remote pipeline to a specific version, branch, or commit for reproducible execution, and that the pin applies only to that one invocation, not to later runs.

### What's next?

You've covered the fundamentals of running and managing Nextflow pipelines.
See [Course summary](next_steps.md) for where to go from here.

---

## Summary

In this part you learned to:

- Run a pipeline directly from a GitHub repository without downloading it
- Pin a remote pipeline to a specific revision for reproducibility, and understand that the pin applies only to that one invocation
