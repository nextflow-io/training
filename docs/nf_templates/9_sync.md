# TEMPLATE syncs

The template evolves as the ecosystem evolves.

New templates are release semi-regularly and the nf-core tooling helps incorporate these changes through syncs with the `TEMPLATE` branch.

```bash
git add .
git commit -m "Added fastp"
git push
```

```bash
nf-core sync
```

The tooling merges updates suggesting a git command.

```bash
cd /Users/chris/workspace/myorg-mypipeline                            
git merge TEMPLATE 
```

## Bump the version across

The pipeline version can be bumped across the template using nf-core tooling.

```bash
nf-core bump-version 1.0
```

## Update GitHub

Changes should be pushed to GitHub.

```bash
git add .
git commit -m "Version release"
git push
```
