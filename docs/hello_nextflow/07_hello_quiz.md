# Hello Nextflow Quiz

Great stuff, you've completed all of the training material in Hello Nextflow!
But, did it stick?

Test yourself with the quiz below to see how many Nextflow fundamentals you can remember :think:

<?quiz?>

question: Are you ready?
answer-correct: Yes!
answer: No!
answer: Maybe!
content:

<?/quiz?>

## Part 1: Hello World

<?quiz?>

question: What is the minimum requirement for a Nextflow process?
answer: Input and output blocks
answer-correct: Output and script blocks
answer: Only a script block
answer: Input block
content:

<?/quiz?>

<?quiz?>

question: What does the output block in a process do?
answer: Creates the output files automatically
answer-correct: Declares what output files to expect from the process
answer: Determines the format of the output
answer: Publishes outputs to a directory
content:

<?/quiz?>

<?quiz?>

question: What command is used to run a Nextflow workflow?
answer: `nextflow execute script.nf`
answer-correct: `nextflow run script.nf`
answer: `nextflow start script.nf`
answer: `nextflow launch script.nf`
content:

<?/quiz?>

<?quiz?>

question: What is stored in the `work/` directory?
answer: Only the final output files
answer: Configuration files
answer-correct: Temporary files for each process execution
answer-correct: Execution metadata like `.command.sh` and `.exitcode`
answer-correct: Output files from processes
answer: Container images
content:

<?/quiz?>

<?quiz?>

question: What does the `-resume` flag do?
answer: Restarts the entire workflow from scratch
answer-correct: Skips processes that already ran with the same code, settings, and inputs
answer: Resumes only failed processes
answer: Saves the workflow state for later
content:

<?/quiz?>

<?quiz?>

question: What is the default mode for the `publishDir` directive?
answer: copy
answer-correct: symlink
answer: move
answer: hardlink
content:

<?/quiz?>

<?quiz?>

question: How do you pass a parameter to a Nextflow workflow from the command line?
answer: `-parameter value`
answer-correct: `--parameter value`
answer: `params.parameter = value`
answer: `-p parameter=value`
content:

<?/quiz?>

<?quiz?>

question: Which quote type is required for variable interpolation in Nextflow?
answer-correct: Double quotes `"`
answer: Single quotes `'`
answer: Backticks `` ` ``
answer: Triple quotes `"""`
content:

<?/quiz?>

## Part 2: Hello Channels

<?quiz?>

question: What is a channel in Nextflow?
answer: A file directory
answer-correct: A queue-based mechanism for passing data between processes
answer: A configuration parameter
answer: A container registry
content:

<?/quiz?>

<?quiz?>

question: Which channel factory creates a simple queue channel with values?
answer: `channel.fromPath()`
answer-correct: `channel.of()`
answer: `channel.from()`
answer: `channel.create()`
content:

<?/quiz?>

<?quiz?>

question: What happens when you pass a channel with multiple values to a process?
answer: Only the first value is processed
answer-correct: The process is called once for each value in the channel
answer: All values are processed together in a single call
answer: The workflow fails with an error
content:

<?/quiz?>

<?quiz?>

question: What does the `flatten()` operator do?
answer-correct: Unpacks array contents into individual items
answer: Removes duplicate items from a channel
answer: Sorts channel contents alphabetically
answer: Converts items to lowercase
content:

<?/quiz?>

<?quiz?>

question: What is the purpose of the `view()` operator?
answer: To visualize the workflow structure
answer-correct: To inspect and print channel contents for debugging
answer: To create graphical output
answer: To validate channel data types
content:

<?/quiz?>

<?quiz?>

question: What does the `splitCsv()` operator do?
answer: Creates a CSV file from channel data
answer-correct: Parses CSV-formatted content, reading each line into an array
answer: Splits a channel into multiple channels
answer: Validates CSV formatting
content:

<?/quiz?>

<?quiz?>

question: Which operator transforms channel contents element-by-element?
answer: `transform()`
answer: `apply()`
answer-correct: `map()`
answer: `modify()`
content:

<?/quiz?>

<?quiz?>

question: Why do you need dynamic filenames when processing multiple inputs?
answer: To save disk space
answer: To improve performance
answer-correct: To prevent output files from overwriting each other
answer: To make resume work correctly
content:

<?/quiz?>

## Part 3: Hello Workflow

<?quiz?>

question: How do you access the output of a process in Nextflow?
answer: `processName.output`
answer-correct: `processName.out`
answer: `processName.result`
answer: `output.processName`
content:

<?/quiz?>

<?quiz?>

question: What determines the execution order of processes in a workflow?
answer: The order they are defined in the file
answer: Alphabetical order by process name
answer-correct: The dependencies between processes (output → input)
answer: The order they are called in the workflow block
content:

<?/quiz?>

<?quiz?>

question: What does the `collect()` operator do?
answer: Unpacks arrays into individual items
answer-correct: Packages many channel items into a single array element
answer: Removes null values from a channel
answer: Collects error messages
content:

<?/quiz?>

<?quiz?>

question: When should you use the `collect()` operator?
answer: Before processing items individually in parallel
answer-correct: Before a process that needs all results at once
answer: To speed up channel operations
answer: When publishing outputs
content:

<?/quiz?>

<?quiz?>

question: How do you access named outputs from a process?
answer: `processName.out[outputname]`
answer-correct: `processName.out.outputname`
answer: `processName.outputname`
answer: `output.processName.outputname`
content:

<?/quiz?>

<?quiz?>

question: What is the syntax for naming a process output?
answer-correct: `output: path "file.txt", emit: myoutput`
answer: `output: path "file.txt", name: myoutput`
answer: `output: path "file.txt" as myoutput`
answer: `output: path "file.txt" -> myoutput`
content:

<?/quiz?>

<?quiz?>

question: When passing multiple inputs to a process, what must match?
answer-correct: The order of inputs in the workflow call must match the input block order
answer: The data types must all be the same
answer: The variable names must match exactly
answer: The number of inputs must equal the number of outputs
content:

<?/quiz?>

## Part 4: Hello Modules

<?quiz?>

question: What is a Nextflow module?
answer: A configuration file
answer: A workflow definition
answer-correct: A standalone file containing a single process definition
answer: A container image
content:

<?/quiz?>

<?quiz?>

question: What is the recommended filename convention for modules?
answer: `module_processName.nf`
answer-correct: `processName.nf`
answer: `processName.module`
answer: `processName_module.nf`
content:

<?/quiz?>

<?quiz?>

question: Where should module files typically be stored?
answer: In the root directory
answer: In the `lib/` directory
answer-correct: In the `modules/` directory
answer: In the `processes/` directory
content:

<?/quiz?>

<?quiz?>

question: What is the correct syntax to import a module?
answer: `import { processName } from './modules/processName.nf'`
answer-correct: `include { processName } from './modules/processName.nf'`
answer: `require { processName } from './modules/processName.nf'`
answer: `use { processName } from './modules/processName.nf'`
content:

<?/quiz?>

<?quiz?>

question: What happens to resume functionality when you convert processes to modules?
answer: Resume stops working
answer: Resume only works for the main workflow
answer-correct: Resume continues to work normally
answer: Resume works but cache is reset
content:

<?/quiz?>

<?quiz?>

question: What are key benefits of using modules?
answer: Faster execution time
answer: Reduced memory usage
answer-correct: Reuse processes in multiple workflows without duplication
answer-correct: Easier maintenance - improve module once, all workflows benefit
answer-correct: Better code organization
answer: Automatic parallelization
content:

<?/quiz?>

## Part 5: Hello Containers

<?quiz?>

question: What is a container?
answer: A directory structure for organizing code
answer-correct: A lightweight, standalone executable unit with application and dependencies
answer: A Nextflow configuration file
answer: A type of channel
content:

<?/quiz?>

<?quiz?>

question: What is the difference between a container image and instance?
answer: They are the same thing
answer-correct: Image is a template; instance is a running container created from the image
answer: Image is local; instance is remote
answer: Instance is smaller than image
content:

<?/quiz?>

<?quiz?>

question: What does the `-v` flag do in Docker?
answer: Enables verbose output
answer: Sets the Docker version
answer-correct: Mounts a volume to share directories between host and container
answer: Validates the container image
content:

<?/quiz?>

<?quiz?>

question: Why do you need to mount volumes when using containers?
answer: To improve performance
answer: To reduce disk space usage
answer-correct: Because containers are isolated from the host filesystem by default
answer: To enable networking
content:

<?/quiz?>

<?quiz?>

question: How do you specify a container for a Nextflow process?
answer: `docker 'image_uri'` directive
answer-correct: `container 'image_uri'` directive
answer: `image 'image_uri'` directive
answer: `containerImage 'image_uri'` directive
content:

<?/quiz?>

<?quiz?>

question: How do you enable Docker in a Nextflow workflow?
answer: Use the `-with-docker` flag only
answer-correct: Set `docker.enabled = true` in the config file
answer: Docker is enabled by default
answer: Use the `--docker` parameter
content:

<?/quiz?>

<?quiz?>

question: What does Nextflow automatically do when using containers?
answer: Builds the container from a Dockerfile
answer-correct: Pulls the image if not already local
answer-correct: Mounts the work directory into the container
answer-correct: Runs the process script inside the container
answer-correct: Handles cleanup after execution
answer: Installs Docker on the system
answer: Creates a new container image
content:

<?/quiz?>

## Part 6: Hello Config

<?quiz?>

question: What file does Nextflow automatically load for configuration?
answer: `config.nf`
answer: `workflow.config`
answer-correct: `nextflow.config`
answer: `settings.config`
content:

<?/quiz?>

<?quiz?>

question: What is the configuration precedence in Nextflow?
answer: Config file overrides CLI arguments
answer-correct: CLI arguments override config file values
answer: Default values override everything
answer: Profile settings always take precedence
content:

<?/quiz?>

<?quiz?>

question: Can you enable both Docker and Conda in the same configuration?
answer: No, only one can be enabled at a time
answer-correct: Yes, both can be enabled simultaneously
answer: Only if using different profiles
answer: Only on certain operating systems
content:

<?/quiz?>

<?quiz?>

question: What happens when both Docker and Conda are enabled for a process that has both directives?
answer: Conda is prioritized
answer-correct: Docker is prioritized
answer: The workflow fails with an error
answer: The first one defined is used
content:

<?/quiz?>

<?quiz?>

question: What is the default memory allocation for a Nextflow process?
answer: 1.GB
answer-correct: 2.GB
answer: 4.GB
answer: It depends on the system
content:

<?/quiz?>

<?quiz?>

question: How do you set resources for a specific process in the config?
answer-correct: `process { withName: 'processName' { cpus = 4 } }`
answer: `process.processName.cpus = 4`
answer: `processName { cpus = 4 }`
answer: `resources.processName.cpus = 4`
content:

<?/quiz?>

<?quiz?>

question: What command generates a resource profiling report?
answer: `nextflow run script.nf -profile`
answer: `nextflow run script.nf -report`
answer-correct: `nextflow run script.nf -with-report filename.html`
answer: `nextflow run script.nf --resource-report`
content:

<?/quiz?>

<?quiz?>

question: What is the purpose of `resourceLimits` in the configuration?
answer: To set minimum resource requirements
answer-correct: To cap the maximum resources that can be requested
answer: To automatically allocate optimal resources
answer: To track resource usage over time
content:

<?/quiz?>

<?quiz?>

question: What is the default executor in Nextflow?
answer: `slurm`
answer-correct: `local`
answer: `cloud`
answer: `batch`
content:

<?/quiz?>

<?quiz?>

question: How do you use a parameter file with Nextflow?
answer: `nextflow run script.nf -config file.json`
answer-correct: `nextflow run script.nf -params-file file.json`
answer: `nextflow run script.nf --parameters file.json`
answer: Parameter files are loaded automatically
content:

<?/quiz?>

<?quiz?>

question: What can profiles be used for in Nextflow?
answer: User account settings
answer-correct: Grouping infrastructure settings (executor, container technology)
answer-correct: Setting different resource limits for different environments
answer-correct: Providing test parameter values
answer: Performance benchmarking
answer: Container registry authentication
content:

<?/quiz?>

<?quiz?>

question: How do you use multiple profiles simultaneously?
answer: You can only use one profile at a time
answer: `nextflow run script.nf -profile profile1 -profile profile2`
answer-correct: `nextflow run script.nf -profile profile1,profile2`
answer: Profiles are automatically merged
content:

<?/quiz?>

## General Concepts

<?quiz?>

question: Which of these is part of the Nextflow language?
answer-correct: Process
answer: Rule
answer-correct: Workflow
answer-correct: Module
answer: Switch
answer-correct: Channel
answer-correct: Include
content:

<?/quiz?>

<?quiz?>

question: What does DSL2 refer to?
answer: Docker Specification Language version 2
answer-correct: Domain-Specific Language version 2
answer: Data Structure Library version 2
answer: Distributed System Layer version 2
content:

<?/quiz?>

<?quiz?>

question: What is the shebang required for Nextflow scripts?
answer: `#!/bin/bash`
answer-correct: `#!/usr/bin/env nextflow`
answer: `#!/usr/bin/nextflow`
answer: No shebang is required
content:

<?/quiz?>

<?quiz?>

question: What is the naming convention for process names in DSL2?
answer: lowercase
answer: camelCase
answer-correct: UPPERCASE
answer: snake_case
content:

<?/quiz?>

<?quiz?>

question: Which tool can generate container images from Conda packages?
answer: DockerHub
answer-correct: Seqera Containers
answer: BioContainers
answer: Conda Cloud
content:

<?/quiz?>

<?quiz?>

question: What type of quotes should be used for multi-line script blocks?
answer: Single quotes `'...'`
answer: Double quotes `"..."`
answer-correct: Triple quotes `"""..."""`
answer: Backticks `` `...` ``
content:

<?/quiz?>

<?quiz?>

question: What is the file extension for Nextflow scripts?
answer: `.nextflow`
answer-correct: `.nf`
answer: `.nx`
answer: `.workflow`
content:

<?/quiz?>

<?quiz?>

question: Which files in the work directory are useful for debugging? (Select all that apply)
answer-correct: `.command.sh` - shows the actual command executed
answer-correct: `.command.err` - shows stderr output
answer-correct: `.command.out` - shows stdout output
answer-correct: `.command.log` - shows complete log
answer-correct: `.exitcode` - shows the exit status
answer: `.command.config` - shows process configuration
content:

<?/quiz?>

<?quiz?>

question: What happens if you use single quotes for a string that needs variable interpolation?
answer: Nextflow automatically converts them to double quotes
answer-correct: The variable name appears literally; interpolation doesn't work
answer: The workflow fails with a syntax error
answer: The variable is interpolated anyway
content:

<?/quiz?>

<?quiz?>

question: What does the `nextflow clean` command do?
answer: Cleans up container images
answer-correct: Removes work subdirectories from past runs
answer: Deletes configuration files
answer: Clears the Nextflow cache
content:

<?/quiz?>
