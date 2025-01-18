#!/bin/bash

# Define variables
prefix="HybridAssembly_Comp50Cont5"
work_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/magscreen
output_dir=${work_dir}/Results/${prefix}

# Get list of novel MAGs
cd ${output_dir}/new_species
ls -1 *fa | sort -V > ../new_species.txt
cd ${work_dir}

# Add taxonomic metadata to novel MAGs/new species identified from magscreen
Rscript 2.Tax_annot_novel_MAGs.R

