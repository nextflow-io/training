
cd /workspace/gitpod/nf-training/variantcalling


### run through with single sample calling

nextflow run nf-core/sarek \
--input /workspace/gitpod/nf-training/data/reads/variantcalling/sarek-input.csv \
--outdir . \
--tools haplotypecaller \
--genome GRCh38chr21 \
--skip_tools haplotypecaller_filter

## run through with joint calling

nextflow run nf-core/sarek \
--input /workspace/gitpod/nf-training/data/reads/variantcalling/sarek-input.csv \
--outdir . \
--tools haplotypecaller \
--genome GRCh38chr21 \
--skip_tools haplotypecaller_filter \
--joint_germline \
--intervals /workspace/gitpod/nf-training/variantcalling/chr21_intervals.list \
-resume

