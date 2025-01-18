#!/bin/bash

# Move and rename QC'ed reads into a new folder
readQC_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/1.READ_QC"
cleaned_reads_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/1.5.CLEAN_READS"

for i in ${readQC_dir}/*; do 
        sample=$(basename "${i}")
        mv ${i}/final_pure_reads_1.fastq ${cleaned_reads_dir}/${sample}_1.fastq
        mv ${i}/final_pure_reads_2.fastq ${cleaned_reads_dir}/${sample}_2.fastq
done
