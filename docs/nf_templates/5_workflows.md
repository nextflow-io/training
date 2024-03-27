# GitHub actions workflows

Automation is an important part of the nf-core pipeline template. They are used to ensure that pipelines are written to a set standard and automate tasks.

Automated tests with GitHub Actions are configured in the `.github/workflows` folder. Each file defines several tests, for example, running the pipeline with the test data, branch protection, and linting the code and documentation.

- `branch.yml`
    - Sets the branch protection for the nf-core repository
    - _Configured for nf-core repo only_
- `ci.yml`
    - Run small pipeline tests with the small test datasets
    - _Configured for all pipelines_
- `clean-up.yml`
    - Automated testing for stale and closed GitHub issues and PRs in nf-core repo
    - _Configured for nf-core repo only_
- `download_pipeline.yml`
    - Test a pipeline download with `nf-core download`.
    - _Configured for nf-core repo only_
- `fix-linting.yml`
    - Fix linting by adding a comment to a PR
    - _Configured for nf-core repo only_
- `linting_comment.yml`
    - Triggered after the linting action and posts an automated comment to the PR, even if the PR is coming from a fork
    - _Configured for nf-core repo only_
- `linting.yml`
    - Triggered on pushes and PRs to the repository and runs `nf-core lint` and markdown lint tests to ensure that the code meets the nf-core guidelines
    - _Configured for all pipelines_
- `release-announcements.yml`
    - Automatic release toot and tweet announcements for nf-core pipeline releases
    - _Configured for nf-core repo only_

Many of these tests are configured for the nf-core repo. However, they can be modified to better suit your needs or ignored.

!!! warning "Deleting workflows"

    Even though many of these workflows are not relevant for private repos. It is recommended to keep them in place to prevent `nf-core lint` from throwing errors. Extra steps or modifications to these workflows can be added for repositories outside of nf-core.
