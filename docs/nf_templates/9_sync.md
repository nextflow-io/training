# TEMPLATE syncs

The template evolves as the ecosystem evolves.

New templates are released semi-regularly and the nf-core tooling helps incorporate these changes through syncs with the `TEMPLATE` branch.

To keep nf-core pipelines up to date with improvements in the main template, we use a method of synchronisation with the `TEMPLATE` branch.

To sync the template, you first need to commit and push your changes to GitHub.

```bash
git add .
git commit -m "Added fastp"
git push
```

This `nf-core sync` command updates the `TEMPLATE` branch with the latest version of the nf-core template, so that these updates can be synchronised with the pipeline. It is run automatically for all pipelines when ever a new release of nf-core/tools (and the included template) is made.

```bash
nf-core sync
```

The tooling merges updates suggesting a git command

```bash
cd /Users/chris/workspace/nf-core-mypipeline
git merge TEMPLATE
```

!!! note

    As this is a newly created template pipeline the `TEMPLATE` branch doesn't need to be synced.
