# Custom modules

nf-core offers a comprehensive set of modules that have been created and curated by the community. However, as a developer, you may be interested in bespoke pieces of software that are not apart of the nf-core repository or customizing a module that already exists.

## Adding a local module

Local modules

## Patching modules

Although nf-core modules are written to be flexible you may want to modify them to better fit your purpose.

The `nf-core lint` command will help manage nf-core components and test that they match the remote source they came from.

For example, if you modify an nf-core module, it will no longer match the remote and a linting test of this module will fail.

!!! question "Exercise"

    Edit the `SEQTK_TRIM` module by adding an `ed` to the end of `reads_fail`. Check to see if your change has caused the linting test to fail.

    ```groovy title="modules/nf-core/fastp/main.nf" linenums="22"
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_failed
    ```
