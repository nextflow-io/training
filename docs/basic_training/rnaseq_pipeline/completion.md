# Pipeline completion

## Handle completion event

This step shows how to execute an action when the pipeline completes the execution.

Note that Nextflow processes define the execution of **asynchronous** tasks i.e. they are not executed one after another as if they were written in the pipeline script in a common **imperative** programming language.

The script uses the `workflow.onComplete` event handler to print a confirmation message when the script completes.

Try to run it by using the following command:

```bash
nextflow run script7.nf -resume --reads 'data/ggal/*_{1,2}.fq'
```

## Email notifications

Send a notification email when the workflow execution completes using the `-N <email address>` command-line option.

Note: this requires the configuration of a SMTP server in the nextflow config file. Below is an example `nextflow.config` file showing the settings you would have to configure:

```groovy
mail {
    from = 'info@nextflow.io'
    smtp.host = 'email-smtp.eu-west-1.amazonaws.com'
    smtp.port = 587
    smtp.user = "xxxxx"
    smtp.password = "yyyyy"
    smtp.auth = true
    smtp.starttls.enable = true
    smtp.starttls.required = true
}
```

See [mail documentation](https://www.nextflow.io/docs/latest/mail.html#mail-configuration) for details.

## Metrics and reports

Nextflow can produce multiple reports and charts providing several runtime metrics and execution information.

Run the [rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf) pipeline previously introduced as shown below:

```bash
nextflow run rnaseq-nf -with-docker -with-report -with-trace -with-timeline -with-dag dag.png
```

The `-with-docker` option launches each task of the execution as a Docker container run command.

The `-with-report` option enables the creation of the workflow execution report. Open the file `report.html` with a browser to see the report created with the above command.

The `-with-trace` option enables the creation of a tab separated file containing runtime information for each executed task. Check the `trace.txt` for an example.

The `-with-timeline` option enables the creation of the workflow timeline report showing how processes were executed over time. This may be useful to identify the most time consuming tasks and bottlenecks. See an example at [this link](https://www.nextflow.io/docs/latest/tracing.html#timeline-report).

Finally, the `-with-dag` option enables the rendering of the workflow execution direct acyclic graph representation. Note: This feature requires the installation of [Graphviz](http://www.graphviz.org/) on your computer. See [here](https://www.nextflow.io/docs/latest/tracing.html#dag-visualisation) for further details. Then try running :

```bash
dot -Tpng dag.dot > graph.png
open graph.png
```

!!! warning

    Run time metrics may be incomplete for run short running tasks as in the case of this tutorial.

!!! info

    You view the HTML files by right-clicking on the file name in the left side-bar and choosing the **Preview** menu item.
