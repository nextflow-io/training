process KNIT_PHYLOSEQ {
    tag "knit_phyloseq"
    publishDir "$params.outdir", mode:'copy'
    container "community.wave.seqera.io/library/bioconductor-phyloseq_knit_r-base_r-ggplot2_r-rmdformats:6efceb52eb05eb44"

    input:
    path "merged"

    output:
    stdout

    script:
    def report = params.report
    def outdir = params.outdir
    """
    biom_path=\$(realpath ${merged})
    outreport=\$(realpath ${outdir})
    Rscript -e "rmarkdown::render('${report}', params=list(args='\${biom_path}'),output_file='\${outreport}/report.html')"
    """
}
