#!/bin/bash
# runs mutect2 on the specified chromosome for all pairs of files in provided in file_pairs.txt
# usage: ./call_variants_jung_rnaseq.sh chr1
# outputs: vcf.gz files in output_dir, e.g. SRR1014768_vs_SRR1014769_chr1.vcf.gz

# take in a command line argument of chrom number
chrom_name=$1 # e.g., chr1 or chrY
file_pairs_fn=$2 # e.g., /cellar/users/zkoch/dream/data/jung_2015/sample_name_mapping_cleaned.tsv
pair_num=$3 # between 1 and 10
bam_dir=/cellar/users/zkoch/dream/data/jung_2015/STAR_results_two_pass
output_dir=/cellar/users/zkoch/dream/data/jung_2015/called_variants_bsqrd
# path to reference files
reference_file=/cellar/users/zkoch/dream/data/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa
# these are from tutorial followed in somatic mutation projects /cellar/users/zkoch/methylation_and_mutation/data/normal_try2_tcga/mutation_calling/tutorial_files
# they say chr17 but aren't actually
pon_fn=/cellar/users/zkoch/methylation_and_mutation/data/normal_try2_tcga/mutation_calling/tutorial_files/resources/chr17_m2pon.vcf.gz
#nomad_fn="/cellar/users/zkoch/methylation_and_mutation/data/normal_try2_tcga/mutation_calling/tutorial_files/resources/chr17_af-only-gnomad_grch38.vcf.gz"
gnomad_fn="/cellar/users/zkoch/methylation_and_mutation/data/normal_try2_tcga/mutation_calling/tutorial_files/resources/gnomad_grch38.renamedChrs.vcf.gz"

# read run_accession_early column of file_pairs_fn into a bash array
readarray -t early_names < <(awk -F'\t' '{print $2}' $file_pairs_fn)
readarray -t late_names < <(awk -F'\t' '{print $3}' $file_pairs_fn)

# loop through all pairs
#for ((pair_num=1; pair_num<${#early_names[@]}; pair_num++))
#do
# get the file names that correspond to ${late_names[pair_num]} and ${early_names[pair_num]}
# early_fn=$bam_dir/${early_names[pair_num]}Aligned.sortedByCoord.out.bam.marked_duplicates.split.named.bam
# late_fn=$bam_dir/${late_names[pair_num]}Aligned.sortedByCoord.out.bam.marked_duplicates.split.named.bam

early_fn=$bam_dir/${early_names[pair_num]}Aligned.sortedByCoord.out.bam.marked_duplicates.split.named.recalibrated.bam
late_fn=$bam_dir/${late_names[pair_num]}Aligned.sortedByCoord.out.bam.marked_duplicates.split.named.recalibrated.bam

out_fn="$output_dir/${late_names[pair_num]}_vs_${early_names[pair_num]}_chr${chrom_name}.vcf.gz"
# print these
echo "pair_num": $pair_num
echo "early_fn": $early_fn
echo "late_fn": $late_fn
echo "out_fn": $out_fn
# check if $out_fn exists already
if [ -f $out_fn ]; then
    echo "File" $out_fn "already exists, skipping"
else
    echo "Running Mutect2"
    # run mutec
    gatk Mutect2 \
        -R $reference_file \
        -I $late_fn \
        -I $early_fn \
        -normal ${early_names[pair_num]} \
        -O $out_fn \
        --germline-resource $gnomad_fn \
        -L $chrom_name
        --pcr-snv-qual 20 \
        --pcr-indel-qual 20 \
        --callable-depth 4 \
        #-pon $pon_fn \
fi
#done