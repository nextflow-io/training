# Instructor notes: Use nf-core

Notes for instructors and teaching assistants running the [Use nf-core](https://training.nextflow.io/latest/nfcore_use/) course live.
Read this before the session; learners don't need it.

## Before the session

- **Use the right version of the materials.**
  If you're teaching content that is on `master` but not yet released, point learners to the development version of the website (version menu), not `latest`.
  Released versions open Codespaces on the release tag, so their files may not match unreleased lessons.
- **Have learners create a new codespace.**
  The setup link uses `quickstart=1`, which reopens an existing codespace if a learner already has one for this repository, with its old files.
  At the start, either give out a link without it (`https://codespaces.new/nextflow-io/training?ref=master`), or tell learners to create a new codespace if GitHub offers to resume one.
  If someone is already in an old codespace, `git pull` in `/workspaces/training` brings it up to date.
- **Check the environment** once the codespace is up: `ls nfcore-use` should show `data/`, `broken_samplesheet_1.csv` and `solutions/`, and `nextflow -version` should match the version in the orientation page.
- **Don't merge changes to the course branch during the session**, so the site and codespaces stay the same as what you tested.

## Machine resources

- The default codespace has **2 CPUs, about 8 GB of memory and a 32 GB disk**.
  Parts 1 to 5 run on it.
- Part 6 asks learners to switch to the **8-core machine (32 GB memory, 64 GB disk)**.
  This uses their quota 4x faster (about 15 hours of the free 120 core-hours per month), restarts the codespace (a few minutes, files are kept, running commands stop), and learners need to `cd /workspaces/training/nfcore-use` again afterwards.
  Remind learners to switch back at the end.
- **Disk space:** a full `nf-core/rnaseq` test run on the default machine can fill the disk with container images and fail with `no space left on device`.
  The course has learners cancel that run early, but if someone lets it finish: `docker image prune -a`.

## Expected errors and confusing output, by part

### Parts 1 and 2

- Values passed on the command line arrive as **strings** in recent Nextflow versions, for booleans and integers alike, so `--skip_trim true` fails schema validation.
  The course uses params files for these; that is the fix if a learner tries the command line.
- The help message suggests `--showHidden`, but only `--show_hidden` works.

### Part 3

- The first `nf-core/rnaseq` run (section 1.2) **fails on purpose** with `Process requirement exceeds available memory -- req: 12 GB`.
- `laptop.config` caps CPUs as well as memory.
  `TRIMGALORE` sets its CPUs by full process name in the pipeline, so it needs its own `withName: '.*:TRIMGALORE'` override; a selector using only `TRIMGALORE` does not take effect.
  If a learner sees `Process requirement exceeds available CPUs -- req: 4; avail: 2`, they are using an old `laptop.config`.
- After Ctrl+C, nf-core prints `Pipeline completed with errors` and `Pipeline failed`.
  This is expected when a run is cancelled.
  Learners should press Ctrl+C once and wait a few seconds.

### Part 4

- Nothing runs in the codespace; everything happens in the browser on nf-co.re.
- Each scenario has a collapsed hint listing candidate pipelines.
  The lesson suggests opening it after about five minutes; in a group session, about 10 minutes works better, since groups spend time discussing.
- Some published test results are incomplete: recent `nf-core/raredisease` releases only publish `pipeline_info/`, and the `nf-core/longraredisease` test covers one chromosome of one sample, without the AnnotSV or SvAnna outputs.
  Point learners to the output documentation for those.
- Worked answers are in the collapsed solution blocks at the end of the page.

### Part 5

- `nf-core/fetchngs` downloads from the ENA, so it needs internet access to `ftp.sra.ebi.ac.uk`.
  If the network blocks it, demo the step from your screen; learners can follow section 3.2 from the generated samplesheet shown in the lesson.
- The fetchngs and `-stub` commands include `-c laptop.config`: without it, they fail on the default machine with the same 12 GB memory error as in Part 3.
- `-preview` only triggers schema validation; the cross-row check (mismatched strandedness) only appears with `-stub`.
- When the fixed samplesheet passes the `-stub` run, tasks start running; learners can stop it with Ctrl+C.

### Part 6

- The default-requests run takes about 14 minutes on a new 8-core codespace, mostly downloading container images.
  In a timed session, have learners start it before a break.
- The run with `tuned.config` takes about 5.5 minutes, compared to about 6.5 minutes with the default requests once the images are downloaded.
  The time saved on this tiny test dataset is modest; the point is the memory reservation and the number of tasks running in parallel, visible in the execution timeline.
- Learners' resource numbers will differ slightly from the tables in the lesson.

## Suggested timing

| Part                                   | Duration                                                 |
| -------------------------------------- | -------------------------------------------------------- |
| 1. Run a demo pipeline                 | 20 min                                                   |
| 2. Configure pipeline execution        | 20 min                                                   |
| 3. Run a production pipeline           | 25 min                                                   |
| 4. Choose a pipeline                   | 40 min                                                   |
| 5. Prepare your inputs                 | 45 min                                                   |
| 6. Scale up and tune compute resources | 40 min, including a machine switch and two pipeline runs |

## Known upstream issues

These may come up if learners explore further, or ask why the course works around them:

- `nf-core/fetchngs` 1.13.0 requests 12 GB and 2 CPUs per download (up to 36 GB on retry) for a task that uses about 15 MB; this is why Part 5 runs it with `-c laptop.config`.
- `nf-core/differentialabundance` 2.0.0: the Quarto report step fails on Seqera Platform with Fusion enabled (`attempt to write a readonly database`); the analysis itself completes.
- `nf-core/raredisease` 3.1.2: `PREP_MITOSALT` has no container and fails on AWS Batch.
