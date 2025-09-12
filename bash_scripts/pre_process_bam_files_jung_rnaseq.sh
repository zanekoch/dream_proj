#!/bin/bash
# runs mutect2 on the specified chromosome for all pairs of files in provided in file_pairs.txt
# usage: ./call_variants_jung_rnaseq.sh chr1
# outputs: vcf.gz files in output_dir, e.g. SRR1014768_vs_SRR1014769_chr1.vcf.gz

file_pairs_fn=$1 # e.g., /cellar/users/zkoch/dream/data/jung_2015/sample_name_mapping_cleaned.tsv
pair_num=$2 # between 1 and 10
bam_dir=/cellar/users/zkoch/dream/data/jung_2015/STAR_results_two_pass
output_dir=/cellar/users/zkoch/dream/data/jung_2015/called_variants
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
early_fn=$bam_dir/${early_names[pair_num]}Aligned.sortedByCoord.out.bam
late_fn=$bam_dir/${late_names[pair_num]}Aligned.sortedByCoord.out.bam
# print these
echo "pair_num": $pair_num
echo "early_fn": $early_fn
echo "late_fn": $late_fn
# check if $early_fn.marked_duplicates.bam and $late_fn.marked_duplicates.bam exist already
if [ -f $early_fn.marked_duplicates.bam ] && [ -f $late_fn.marked_duplicates.bam ]; then
    echo "Duplicate Files already exist, skipping"
    echo $early_fn.marked_duplicates.bam
    echo $late_fn.marked_duplicates.bam
else
    echo "Running MarkDuplicates"
    # mark duplicates
    gatk MarkDuplicates \
        -I $early_fn \
        -O $early_fn.marked_duplicates.bam \
        -M $early_fn.marked_dup_metrics.txt 
    gatk MarkDuplicates \
        -I $late_fn \
        -O $late_fn.marked_duplicates.bam \
        -M $late_fn.marked_dup_metrics.txt
fi
# add name to .marked_duplicates.bam files
# check if named files exist already
if [ -f $early_fn.marked_duplicates.named.bam ] && [ -f $late_fn.marked_duplicates.named.bam ]; then
    echo "Named Files already exist, skipping"
    echo $early_fn.marked_duplicates.named.bam
    echo $late_fn.marked_duplicates.named.bam
else
    echo "Running AddOrReplaceReadGroups"
    # add read groups, for the purpose of naming
    gatk AddOrReplaceReadGroups \
        -I $early_fn.marked_duplicates.bam \
        -O $early_fn.marked_duplicates.named.bam \
        -LB ${early_names[pair_num]} \
        -PL illumina \
        -PU ${early_names[pair_num]} \
        -SM ${early_names[pair_num]}
    gatk AddOrReplaceReadGroups \
        -I $late_fn.marked_duplicates.bam \
        -O $late_fn.marked_duplicates.named.bam \
        -LB ${late_names[pair_num]} \
        -PL illumina \
        -PU ${late_names[pair_num]} \
        -SM ${late_names[pair_num]}
    # index these files
    samtools index $early_fn.marked_duplicates.named.bam $early_fn.marked_duplicates.named.bam.bai
    samtools index $late_fn.marked_duplicates.named.bam $late_fn.marked_duplicates.named.bam.bai        
fi
# recalibrate base quality scores
# check if recalibrated files exist already
if [ -f $early_fn.marked_duplicates.named.recalibrated.bam ] && [ -f $late_fn.marked_duplicates.named.recalibrated.bam ]; then
    echo "Recalibrated Files already exist, skipping"
    echo $early_fn.marked_duplicates.named.recalibrated.bam
    echo $late_fn.marked_duplicates.named.recalibrated.bam
else
    echo "Running BaseRecalibrator"
    # base recalibrator
    gatk BaseRecalibrator \
        -R $reference_file \
        -I $early_fn.marked_duplicates.named.bam \
        --known-sites $gnomad_fn \
        -O $early_fn.marked_duplicates.named.recal_data.table
    gatk BaseRecalibrator \
        -R $reference_file \
        -I $late_fn.marked_duplicates.named.bam \
        --known-sites $gnomad_fn \
        -O $late_fn.marked_duplicates.named.recal_data.table
    echo "Running ApplyBQSR"
    # apply bqsr
    gatk ApplyBQSR \
        -R $reference_file \
        -I $early_fn.marked_duplicates.named.bam \
        --bqsr-recal-file $early_fn.marked_duplicates.named.recal_data.table \
        -O $early_fn.marked_duplicates.named.recalibrated.bam
    gatk ApplyBQSR \
        -R $reference_file \
        -I $late_fn.marked_duplicates.named.bam \
        --bqsr-recal-file $late_fn.marked_duplicates.named.recal_data.table \
        -O $late_fn.marked_duplicates.named.recalibrated.bam
    # index these files
    samtools index $early_fn.marked_duplicates.named.recalibrated.bam $early_fn.marked_duplicates.named.recalibrated.bam.bai
    samtools index $late_fn.marked_duplicates.named.recalibrated.bam $late_fn.marked_duplicates.named.recalibrated.bam.bai
fi

# check if $early_fn.marked_duplicates.split.bam and $late_fn.marked_duplicates.split.bam exist already
if [ -f $early_fn.marked_duplicates.split.bam ] && [ -f $late_fn.marked_duplicates.split.bam ]; then
    echo "Split Files already exist, skipping"
    echo $early_fn.marked_duplicates.split.bam
    echo $late_fn.marked_duplicates.split.bam
else
    echo "Running SplitNCigarReads"
    # splitNcigarReads, because RNA-seq
    gatk SplitNCigarReads \
        -R $reference_file \
        -I $early_fn.marked_duplicates.bam \
        -O $early_fn.marked_duplicates.split.bam
    gatk SplitNCigarReads \
        -R $reference_file \
        -I $late_fn.marked_duplicates.bam \
        -O $late_fn.marked_duplicates.split.bam
fi
# add names to .bam files
# check if named files exist already
if [ -f $early_fn.marked_duplicates.split.named.bam ] && [ -f $late_fn.marked_duplicates.split.named.bam ]; then
    echo "Named Files already exist, skipping"
    echo $early_fn.marked_duplicates.split.named.bam
    echo $late_fn.marked_duplicates.split.named.bam
else
    echo "Running AddOrReplaceReadGroups"
    # add read groups, for the purpose of naming
    gatk AddOrReplaceReadGroups \
        -I $early_fn.marked_duplicates.split.bam \
        -O $early_fn.marked_duplicates.split.named.bam \
        -LB ${early_names[pair_num]} \
        -PL illumina \
        -PU ${early_names[pair_num]} \
        -SM ${early_names[pair_num]}
    gatk AddOrReplaceReadGroups \
        -I $late_fn.marked_duplicates.split.bam \
        -O $late_fn.marked_duplicates.split.named.bam \
        -LB ${late_names[pair_num]} \
        -PL illumina \
        -PU ${late_names[pair_num]} \
        -SM ${late_names[pair_num]}
    # index these files
    samtools index $early_fn.marked_duplicates.split.named.bam $early_fn.marked_duplicates.split.named.bam.bai
    samtools index $late_fn.marked_duplicates.split.named.bam $late_fn.marked_duplicates.split.named.bam.bai        
fi
# recalibrate the split files also
# check if recalibrated files exist already
if [ -f $early_fn.marked_duplicates.split.named.recalibrated.bam ] && [ -f $late_fn.marked_duplicates.split.named.recalibrated.bam ]; then
    echo "Recalibrated Files already exist, skipping"
    echo $early_fn.marked_duplicates.split.named.recalibrated.bam
    echo $late_fn.marked_duplicates.split.named.recalibrated.bam
else
    echo "Running BaseRecalibrator"
    # base recalibrator
    gatk BaseRecalibrator \
        -R $reference_file \
        -I $early_fn.marked_duplicates.split.named.bam \
        --known-sites $gnomad_fn \
        -O $early_fn.marked_duplicates.split.named.recal_data.table
    gatk BaseRecalibrator \
        -R $reference_file \
        -I $late_fn.marked_duplicates.split.named.bam \
        --known-sites $gnomad_fn \
        -O $late_fn.marked_duplicates.split.named.recal_data.table
    echo "Running ApplyBQSR"
    # apply bqsr
    gatk ApplyBQSR \
        -R $reference_file \
        -I $early_fn.marked_duplicates.split.named.bam \
        --bqsr-recal-file $early_fn.marked_duplicates.split.named.recal_data.table \
        -O $early_fn.marked_duplicates.split.named.recalibrated.bam
    gatk ApplyBQSR \
        -R $reference_file \
        -I $late_fn.marked_duplicates.split.named.bam \
        --bqsr-recal-file $late_fn.marked_duplicates.split.named.recal_data.table \
        -O $late_fn.marked_duplicates.split.named.recalibrated.bam
    # index these files
    samtools index $early_fn.marked_duplicates.split.named.recalibrated.bam $early_fn.marked_duplicates.split.named.recalibrated.bam.bai
    samtools index $late_fn.marked_duplicates.split.named.recalibrated.bam $late_fn.marked_duplicates.split.named.recalibrated.bam.bai
fi