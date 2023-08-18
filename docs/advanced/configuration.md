# Configuration

To be introduced on day #2.

This is an aspect of Nextflow that can be confusing. There are multiple ways of loading configuration and parameters into a Nextflow.

This gives us two complications:

-   At which location should I be loading a configuration value?
-   Given a particular parameter, how do I know where it was set?

## Precedence

1. Parameters specified on the command line (`--something value`)
2. Parameters provided using the `-params-file` option
3. Config file specified using the `-c my_config option`
4. The config file named `nextflow.config` in the current directory
5. The config file named `nextflow.config` in the workflow project directory
6. The config file `$HOME/.nextflow/config`
7. Values defined within the pipeline script itself (e.g. `main.nf`)
