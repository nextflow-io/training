# Next Steps

Congratulations on completing the **Troubleshooting Workflows** mini-course.

---

## 1. Apply the toolkit to your own pipelines

The fastest way to internalise these techniques is to reach for them on your own work.
Next time a pipeline of yours fails, before reaching for the search engine:

- Read the error message all the way through, then find the work directory it cites.
- Walk through `.command.sh`, `.command.err`, `.command.out`, `.exitcode`, and `ls` before changing any code.
- If the error is opaque or absent, run with `-preview` and then `debug true` to confirm what your code actually sees.

You'll find that the toolkit answers most questions on its own.

---

## 2. Build debugging into your workflow design

A few habits make pipelines easier to debug from the start:

- Ship every process with a `stub:` directive so `-stub-run` always works.
- Add a `debug` profile to your `nextflow.config` (real-time output, `cleanup = false`, `maxForks = 1`) so a single flag gets you full visibility.
- Use `.view()` liberally during development. Remove it once channels are stable.

---

## 3. Continue your Nextflow training

If you haven't already, check out our other training material:

- **[Hello Nextflow](../../hello_nextflow/index.md)** — Foundational Nextflow concepts
- **[Hello nf-core](../../hello_nf-core/index.md)** — nf-core pipelines and best practices
- **[nf-test](../nf_test/index.md)** — Add tests so future bugs are caught before they hit production
- **[Other Side Quests](../index.md)** — Deep dives into specific topics

---

## Additional resources

- [Nextflow troubleshooting guide](https://www.nextflow.io/docs/latest/troubleshooting.html) — official troubleshooting documentation
- [Nextflow channels reference](https://www.nextflow.io/docs/latest/channel.html) — channel types and behaviour
- [Process directives reference](https://www.nextflow.io/docs/latest/process.html#directives) — all available process configuration options
- [nf-test](https://www.nf-test.com/) — testing framework for Nextflow pipelines
- [Nextflow Slack community](https://www.nextflow.io/slack-invite.html) — help from the community
- [Seqera Platform](https://seqera.io/platform/) — monitoring and debugging at scale for production workflows
