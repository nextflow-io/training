# nf-core components

The nf-core ecosystem includes a repository of components. These components are written and maintained by the nf-core community. Components are written for diverse users and may require additional configuration to suit different use cases. There are two types of components, modules and subworkflows:

- **Modules** are wrappers around single process that contain tests, environments, and meta data.
- **Subworkflows** are two or more modules that are packaged together and also come with tests and metadata.

Currently, there are more than [1200 modules](https://nf-co.re/modules) and [60 subworkflows](https://nf-co.re/subworkflows) (April 2024) available through nf-core.

Both modules and subworkflows can be listed, installed, updated, removed, and patched using nf-core tooling.

## Modules

Modules came packaged with everything needed to run and test the contained process locally and on the cloud:

```
modules/nf-core/fastqc
├── environment.yml
├── main.nf
├── meta.yml
└── tests
    ├── main.nf.test
    ├── main.nf.test.snap
    └── tags.yml
```

### Installing a module

nf-core tools can be used to install components from the nf-core repository:

```
nf-core modules install
```

You can follow the prompts to find and install the module you are interested in:

```console
? Tool name: fastp
```

Once selected, the tooling will install the module in `modules/nf-core/` folder and give you a suggested line of code for your main workflow file (e.g., `workflows/mypipeline.nf`).

```console
INFO     Installing 'fastp'                                                                                                                              
INFO     Use the following statement to include this module:                                                                                             

include { FASTP } from '../modules/nf-core/fastp/main'  
``` 

!!! question "Exercise"

    Install the `fastp` module to your pipeline.

!!! note "What is fastp?

    Fastp is a...

### Adding a module to your pipeline

Although the modules has been installed in your repository, it is not yet added to your pipeline. The suggested `include` statement needs to be added to your `workflows/mypipeline.nf` file.

```console
include { FASTP } from '../modules/nf-core/fastp/main'
```

!!! question "Exercise"

    Add the `include { FASTP } from '../modules/nf-core/fastp/main'` statement to your `workflows/mypipeline.nf` file.

To implement the module in your pipeline you will need to check what inputs are required. You can view the inputs channels for the module by opening the `/modules/nf-core/fastp/main` file.

```groovy title="/modules/nf-core/fastp/main"
input:
tuple val(meta), path(reads)
path  adapter_fasta
val   save_trimmed_fail
val   save_merged
```

Each module has a meta file which describes the required inputs. The meta file is rendered on the [nf-core website](https://nf-co.re/modules/fastp), or can be viewed using the `nf-core modules info command`.

Code for `workflows/mypipeline.nf`

!!! question "Exercise"

    Add the following code to your `workflows/mypipeline.nf` file to include the fastp process in your pipeline.

```groovy
//
// MODULE: Run Fastp
//
FASTP (
    ch_samplesheet,
    [],
    false,
    false
)
```

!!! note "What is `[]`?"

    Every channel needs to be filled for a process to run. In the example above `[]` was used to fill the channel with an empty list.

Importantly, nf-core components are tracked in `modules.json`. Tracking components improves reproducibility and reporting by giving some certainty about the version of the tool you were using. The `modules.json` will be updated with a hash for the module.

!!! tip "Private modules repositories"

    Private repositories of modules can be accessed using `--git-remote` in `nf-core modules install` command.

### Linting modules

The `nf-core lint` command will help manage nf-core components and test that they match the remote they came from. Although modules are written to be flexible, you may want to modify a component to fit your purpose.

For example, if you modify an nf-core module, it will no longer match the remote and a linting test of this module would fail.

```
EXAMPLE
```

!!! question "Exercise"

    Edit the `fastp` module by adding an `s` to the end of `json`. Check to see if your change has caused the linting test to fail.

```
tuple val(meta), path('*.json')           , emit: jsons
```

Running `nf-core lint` after you have modified an nf-core module will cause it to throw an error.

```
╭─ [✗] 1 Module Test Failed ─────────────────────────────────────────────────────────────────────╮
│              ╷                               ╷                                                 │
│ Module name  │ File path                     │ Test message                                    │
│╶─────────────┼───────────────────────────────┼─────────────────────────────────────────────────│
│ fastp        │ modules/nf-core/fastp/main.nf │ Local copy of module does not match remote      │
│              ╵                               ╵                                                 │
╰────────────────────────────────────────────────────────────────────────────────────────────────╯
```

Changing a module does not mean you can't continue to use that module.

The nf-core modules patch command allows you keep using the nf-core component without needing to make it into a local module and curate it yourself. Instead, it creates a path file that will keep track of the changes you made. If you subsequently update the module using the nf-core tooling, the diff file will be retained. If any subsequent changes to the module conflict with your diff file, you will prompted to resolve the conflicts. 

```
nf-core modules patch
```

The prompt can be followed to patch the `fastp` module.

```
? Module name: fastp
```

A patch file is created in the fastp module directory

```
...
INFO     'modules/nf-core/fastp/tests/main.nf.test.snap' is unchanged
INFO     'modules/nf-core/fastp/tests/tags.yml' is unchanged
INFO     'modules/nf-core/fastp/tests/nextflow.config' is unchanged
INFO     'modules/nf-core/fastp/tests/main.nf.test' is unchanged 
INFO     Patch file of 'modules/nf-core/fastp' written to 'modules/nf-core/fastp/fastp.diff' 
```

!!! question "Exercise"

    Patch the `fastp` module to fix the linting error. 