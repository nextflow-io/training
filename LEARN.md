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