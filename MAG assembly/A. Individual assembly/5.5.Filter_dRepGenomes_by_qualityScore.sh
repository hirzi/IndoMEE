#!/bin/bash

# Set working directory
workdir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC

# Generate list of dereplicated genomes
cd /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.4.DREP/dereplicated_genomes
ls -1 *fa | sort > ${workdir}/dereplicated_genomes.txt

# Filter dereplicated genomes by quality score (S = Completion - 5*Contamination >= 50%)
Rscript /rds/project/rds-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/5.5.Filter_dRepGenomes_by_qualityScore.R

# Make new folder containing only genomes that passed quality score filter (if you want to save space, rather than copying you can do symbolic links)
cd ${workdir}
mkdir dereplicated_genomes_filtered
cat dereplicated_genomes_filtered.txt | tail -n +2 | cut -f1 > dereplicated_genomes_filtered_names.txt
rsync -a /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.4.DREP/dereplicated_genomes --files-from=dereplicated_genomes_filtered_names.txt dereplicated_genomes_filtered
