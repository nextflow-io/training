git clone https://github.com/lescai-teaching/datasets_reference_only.git
git clone https://github.com/lescai-teaching/dataset_exercise_resequencing.git

mkdir sarek_run

cd /workspace/gitpod/nf-training/variantcalling/other_data/dataset_exercise_resequencing/raw_data

echo "patient,sample,lane,fastq_1,fastq_2" >>sarek-input.csv
for localfwd in `ls *_1.fq.gz`
do
sample=${localfwd%_1.fq.gz}
localrev=${localfwd%1.fq.gz}2.fq.gz
current=`pwd`
fwd="${current}/${localfwd}"
rev="${current}/${localrev}"
echo "${sample},${sample},L1,${fwd},${rev}" >>sarek-input.csv
done

/workspace/gitpod/nf-training/variantcalling/other_data/dataset_exercise_resequencing/raw_data/sarek-input.csv



nextflow run nf-core/sarek \
--input /workspace/gitpod/nf-training/variantcalling/other_data/dataset_exercise_resequencing/raw_data/sarek-input.csv \
--outdir . \
--tools haplotypecaller \
--genome GRCh38local \
--skip_tools haplotypecaller_filter

## to run joint calling

nextflow run nf-core/sarek \
--input /workspace/gitpod/nf-training/variantcalling/other_data/dataset_exercise_resequencing/raw_data/sarek-input.csv \
--outdir . \
--tools haplotypecaller \
--genome GRCh38local \
--skip_tools haplotypecaller_filter \
--joint_germline \
--intervals chr21_interval.list

