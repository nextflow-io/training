# Bump your pipeline version

The pipeline version number is mentioned in a lot of different places in nf-core pipelines. This command updates the version for you automatically, so that you don't accidentally miss any. It should be used for each pipeline release, and again for the next development version after release.

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

    Create a tagged version release using GitHub

## Update minimum Nextflow version

The template also includes a minimum Nextflow version that is required for the pipeline to to. You can also change the required version of Nextflow using the `nf-core bump-version` command. 

```bash
nf-core bump-version 23.04.0 --nextflow
```
