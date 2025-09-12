download_methylation.txt:
- go to (https://www.encodeproject.org/matrix/?type=Experiment&assay_title=DNAme+array&biosample_ontology.classification=tissue&biosample_ontology.classification=primary+cell) selection methylation array and tissue and primary cell
- click download and select processed files
- xargs -L 1 curl -O -J -L < download_methylation.txt

rna_expression_report_2024_6_9_18h_48m.tsv: all rna seq
- go to (https://www.encodeproject.org/rnaget-report/?type=RNAExpression&file.biosample_ontology.classification=primary+cell&file.assay_title=total+RNA-seq&file.biosample_ontology.classification=tissue&dataset.replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens)
- select tissue and primary cell and human
- donwload all



