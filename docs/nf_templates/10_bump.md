# Bump your pipeline version

The pipeline version number is mentioned in a lot of different places in nf-core pipelines. This  
tool updates the version for you automatically, so that you don't accidentally miss any.  
Should be used for each pipeline release, and again for the next development version after  
release.

```bash
nf-core bump-version 1.0
```

## Update GitHub

After you have updated the version of you pipeline, your changes should be pushed to GitHub.

```bash
git add .
git commit -m "Version 1.0 release"
git push
```

!!! question "Exercise"

    Update your pipeline and push the changes to GitHub

!!! question "Bonus Exercise"

    Create a

## Update

As well as the pipeline version, you can also change the required version of Nextflow.

```bash
nf-core bump-version 23.04.0 --nextflow
```
