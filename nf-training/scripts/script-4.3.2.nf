FROM mambaorg/micromamba
MAINTAINER Your name <your_email>

RUN \
   micromamba install -y -n base -c defaults -c bioconda -c conda-forge \
      salmon=1.5.1 \
      fastqc=0.11.9 \
      multiqc=1.10.1 \
   && micromamba clean -a -y