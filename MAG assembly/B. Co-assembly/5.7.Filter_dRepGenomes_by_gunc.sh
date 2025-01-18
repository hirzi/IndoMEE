#!/bin/bash

# Set working directory
workdir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.5.FILTER_GUNC

# Filter out dereplicated genomes that show chimerism (infered via gunc)
Rscript /rds/project/rds-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/co_assembly/5.7.Filter_dRepGenomes_by_gunc.R

# Make new folder containing only genomes that passed gunc filter (if you want to save space, rather than copying you can do symbolic links)
cd ${workdir}
mkdir dereplicated_genomes_gunc_finalFiltered
cat dereplicated_genomes_gunc_finalFiltered.txt | tail -n +2 | cut -f2 > dereplicated_genomes_gunc_finalFiltered_names.txt
rsync -a /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.4.DREP/dereplicated_genomes --files-from=dereplicated_genomes_gunc_finalFiltered_names.txt dereplicated_genomes_gunc_finalFiltered
