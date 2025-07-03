# Nextflow Command Reference

## Check the Nextflow Version
```sh
nextflow -v
```

# Run a Nextflow Pipeline 
```sh
nextflow run <File.nf>
```
# Resume a Previous Run
```sh
nextflow run hello-world.nf -resume
```
# View Nextflow Logs
```sh
nextflow log 
```

# Dry Run (Preview):
```sh
nextflow clean -before <runnamefromlog> -n
```
-n is for a dry run (shows what would be deleted).

Force Clean (Delete):
```sh
nextflow clean -before <runnamefromlog> -f
```
-f is for force (actually deletes files).

# Passing a variabe 
```sh
nextflow run hello-world.nf --greeting "Hek varlden"
```
---
# More options 
Trace Report
Creates a tabular text file with detailed info about every process:
```sh 
nextflow run pipeline.nf -with-trace
```

Execution Report (HTML)
Generates a report.html file with pipeline stats.
```sh
nextflow run pipeline.nf -with-report
```
Timeline Visualization (HTML)
Generates a timeline.html showing process flow and duration.
```sh
nextflow run pipeline.nf -with-timeline
```
Execution DAG
Creates a flow chart of your pipeline.
```sh
nextflow run pipeline.nf -with-dag flowchart.png
```
