# GitHub actions workflows

Automated workflows are an important part of the nf-core pipeline template.

By default, the template comes with several automated tests that utilize GitHub Actions, each of which are configured in the `.github/workflows` folder. Tests include running the pipeline with the test data, branch protection, and linting the code and documentation. Brief descriptions about each test and file are outlined below:

-   `branch.yml` (_Configured for nf-core repo only_)
    -   Sets the branch protection for the nf-core repository
-   `ci.yml` (_Configured for all pipelines_)
    -   Run small pipeline tests with the small test datasets
-   `clean-up.yml` (_Configured for nf-core repo only_)
    -   Automated testing for stale and closed GitHub issues and PRs in nf-core repo
-   `download_pipeline.yml` (_Configured for nf-core repo only_)
    -   Test a pipeline download with `nf-core download`.
-   `fix-linting.yml` (_Configured for nf-core repo only_)
    -   Fix linting by adding a comment to a PR
-   `linting_comment.yml` (_Configured for nf-core repo only_)
    -   Triggered after the linting action and posts an automated comment to the PR, even if the PR is coming from a fork
-   `linting.yml` (_Configured for all pipelines_)
    -   Triggered on pushes and PRs to the repository and runs `nf-core lint` and markdown lint tests to ensure that the code meets the nf-core guidelines
-   `release-announcements.yml` (_Configured for nf-core repo only_)
    -   Automatic release toot and tweet announcements for nf-core pipeline releases

Notable, many of these tests are only configured for the nf-core repo. However, they can be modified to better suit your needs or ignored if they are superfluous to your requirements.

!!! warning "Deleting workflows"

    Even though many of these workflows are not relevant for private repositories, it is recommended to keep them in place to prevent `nf-core lint` from throwing errors. Extra steps or modifications to these workflows can be added for repositories outside of nf-core.
