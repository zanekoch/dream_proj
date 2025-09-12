#!/bin/bash
# combine the vcfs of each late_vs_early pair using GatherVcfs, then IndexFeatureFile, then MergeMutectStats
# after this step, the vcfs are ready to be filtered using FilterMutectCalls and VariantsToTable

file_pairs_fn=$1 # e.g., /cellar/users/zkoch/dream/data/jung_2015/sample_name_mapping_cleaned.tsv
vcf_dir=$2 # e.g., /cellar/users/zkoch/dream/data/jung_2015/called_variants
reference_file=/cellar/users/zkoch/dream/data/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa

# get file names
readarray -t early_names < <(awk -F'\t' '{print $2}' $file_pairs_fn)
readarray -t late_names < <(awk -F'\t' '{print $3}' $file_pairs_fn)


# iterate across indexes of late_names, skipping the first
for ((i=1; i<${#late_names[@]}; i++))
do
    late_name=${late_names[$i]}
    # 
    echo "late name" $late_name
    # get files mathcing late_name*vcf.gz into a list
    late_vcfs=`ls ${vcf_dir}/${late_name}_vs_*vcf.gz`
    # print the elements of late_vcfs
    echo "${late_vcfs[@]}"
    # check if combined already exists
    if [ -f ${vcf_dir}/${late_name}_all.vcf.gz ]; then
        echo "${vcf_dir}/${late_name}_all.vcf.gz already exists"
        continue
    else
        # turn into list
        late_vcfs=($late_vcfs)
        
        # add each after -I, like -I file1 -I file2 -I file3
        gatk GatherVcfs --OUTPUT ${vcf_dir}/${late_name}_all.vcf.gz --VERBOSITY ERROR --REORDER_INPUT_BY_FIRST_VARIANT -I ${late_vcfs[0]} -I ${late_vcfs[1]} -I ${late_vcfs[2]} -I ${late_vcfs[3]} -I ${late_vcfs[4]} -I ${late_vcfs[5]} -I ${late_vcfs[6]} -I ${late_vcfs[7]} -I ${late_vcfs[8]} -I ${late_vcfs[9]} -I ${late_vcfs[10]} -I ${late_vcfs[11]} -I ${late_vcfs[12]} -I ${late_vcfs[13]} -I ${late_vcfs[14]} -I ${late_vcfs[15]} -I ${late_vcfs[16]} -I ${late_vcfs[17]} -I ${late_vcfs[18]} -I ${late_vcfs[19]} -I ${late_vcfs[20]} -I ${late_vcfs[21]} -I ${late_vcfs[22]} -I ${late_vcfs[23]}
        # create index for late_vcf_fn
        gatk IndexFeatureFile --input ${vcf_dir}/${late_name}_all.vcf.gz
        # create a stats file
        late_stats=`ls ${vcf_dir}/${late_name}_vs_*vcf.gz.stats`
        late_stats=($late_stats)
        gatk MergeMutectStats --stats ${late_stats[0]} --stats ${late_stats[1]} --stats ${late_stats[2]} --stats ${late_stats[3]} --stats ${late_stats[4]} --stats ${late_stats[5]} --stats ${late_stats[6]} --stats ${late_stats[7]} --stats ${late_stats[8]} --stats ${late_stats[9]} --stats ${late_stats[10]} --stats ${late_stats[11]} --stats ${late_stats[12]} --stats ${late_stats[13]} --stats ${late_stats[14]} --stats ${late_stats[15]} --stats ${late_stats[16]} --stats ${late_stats[17]} --stats ${late_stats[18]} --stats ${late_stats[19]} --stats ${late_stats[20]} --stats ${late_stats[21]} --stats ${late_stats[22]} --stats ${late_stats[23]} --output ${vcf_dir}/${late_name}_all.vcf.gz.stats
    fi

    # check if this succeeded
    if [ ! -f "${vcf_dir}/${late_name}_all.vcf.gz" ]; then
        echo "Failed to create ${vcf_dir}/${late_name}_all.vcf.gz"
        exit 1
    fi

    # now filter and create table

    # check if filtered_sample_vcf_fn already exists
    if [ -f "${vcf_dir}/${late_name}_all_filtered.vcf.gz" ]; then
        echo "Filtered VCF file already exists: ${vcf_dir}/${late_name}_all_filtered.vcf.gz, skipping..."
    else
        # filter variants using FilterMutectCalls
        gatk FilterMutectCalls --reference $reference_file --variant ${vcf_dir}/${late_name}_all.vcf.gz --output ${vcf_dir}/${late_name}_all_filtered.vcf.gz
    fi
    # create table of variants that passed the filtering (only variants with PASS in FILTER column)
    if [ -f "${vcf_dir}/${late_name}.tsv" ]; then
        echo "Filtered VCF file already exists: ${vcf_dir}/${late_name}.tsv, skipping..."
    else
        # for VAF calculation
        # DP is read depth
        # AD is alternate counts
        gatk VariantsToTable --reference $reference_file --variant ${vcf_dir}/${late_name}_all_filtered.vcf.gz --output ${vcf_dir}/${late_name}_filtered.tsv -F CHROM -F POS -F REF -F ALT -GF AD -GF DP
        gatk VariantsToTable --reference $reference_file --variant ${vcf_dir}/${late_name}_all.vcf.gz --output ${vcf_dir}/${late_name}_unfiltered.tsv -F CHROM -F POS -F REF -F ALT -GF AD -GF DP

    fi
done