# nf-training-public

Nextflow Training material. 

## Opening the training material in a browser

To open the training, follow the following URL:

https://training.seqera.io/

## Gitpod tutorial

To open the training with a IDE terminal with all programs installed and training material, click the following URL:

https://gitpod.io/#https://github.com/seqeralabs/nf-training-public/tree/master

Then you will need to login to either Github, GitLab or Bitbucket. This will give you 50 hours per month to run Gitpod, with up to 4 parellel workspaces at a time.

Once you have signed in, gitpod should load:

image::gitpod.png[]

Then the IDE should open up and look similar to:

image::gitpod.welcome.png[]

**The sidebar** allows you to customise your environment and perform basic tasks (Copy/Paste, Open files, search, git, etc.) Click the Explorer button to see which files are in this repository:

**The terminal** allows you to run all the programs in the repository, for example `nextflow`, `nf-core` and `docker` are installed in the nf-core rnaseq repository. The terminal may not appear automatically, in which case navigate to the top icon of the sidebar, and choose Terminal/New Terminal.

**The main window** allows you to view and edit files. Clicking on a file in the explorer will open it within the main window. Once a file is open, Markdown or HTML can be rendered using the preview option.



## Rendering ascii docs for training locally

To render the pages on your local machine

First install ascii doctor tools: https://asciidoctor.org/

Then git clone the repo : `git clone --single-branch --branch master https://github.com/seqeralabs/nf-training-public`

Make your changes to the docs in `asciidocs/`.

Then run `make clean`

Then make the docs: `make docs`

Then to open the html: `open asciidocs/index.html`
