#!/bin/bash

work_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metamap
#MAG_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.BIN_REFINEMENT/OLD/Dereplicated_bins_individual_assembly/metawrap_70_5_bins
#MAG_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered
#MAG_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered
#MAG_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.5.FILTER_GUNC_sANI98/dereplicated_genomes_gunc_finalFiltered
#MAG_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP/Comp50Cont5/dereplicated_genomes
#MAG_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP/Comp70Cont5/dereplicated_genomes
#MAG_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP_sfastANI_98/Comp50Cont5/dereplicated_genomes
#MAG_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Public_genDiversity_databases/SPMP/SPMPPolished
MAG_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Public_genDiversity_databases/SPMP/dRep/dereplicated_genomes
#MAG_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_UHGG_ASSEMBLY/5.4.DREP/dereplicated_genomes
#MAG_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_UHGG_SPMP_ASSEMBLY/5.4.DREP/dereplicated_genomes

cd ${MAG_dir}
# Make reference fasta by concatenating all MAGS and renaming constituent contigs based on bin name/number and contig order. 1-indexed.

# Old, derepreciated
#for j in {000000001..000001169}; do 
#	k=$((10#$j))
#	cat bin.$k.fa | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/concat_metawrapIndAssembly_draft1.fa
#done

# New
counter=1
#for j in {000000001..000001162}; do
#for j in {000000001..000000886}; do
#for j in {000000001..000001069}; do
#for j in {000000001..000000988}; do
#for j in {000000001..000004497}; do
#for j in {000000001..000003083}; do
#for j in {000000001..000004629}; do
#for j in {000000001..000004717}; do
#for j in {000000001..000004941}; do
#for j in {000000001..000004860}; do
#for j in {000000001..000001304}; do
#for j in {000000001..000003258}; do
#for j in {000000001..000001210}; do
for j in {000000001..000000646}; do
	#bin=`sed -n ${counter}p < /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered_names.txt`
	#bin=`sed -n ${counter}p < /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered_names.txt`
	#bin=`sed -n ${counter}p < /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.5.FILTER_GUNC_sANI98/dereplicated_genomes_gunc_finalFiltered_names2.txt`
	#bin=`sed -n ${counter}p < /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP/Comp50Cont5/dereplicated_genomes.txt`
	#bin=`sed -n ${counter}p < /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP/Comp70Cont5/dereplicated_genomes.txt`
	#bin=`sed -n ${counter}p < /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP_sfastANI_98/Comp50Cont5/dereplicated_genomes.txt`
	#bin=`sed -n ${counter}p < /home/hl636/rds/rds-mobile-Gwxt74iudAM/Public_genDiversity_databases/SPMP/SPMPPolishedMAGs.txt`
	bin=`sed -n ${counter}p < /home/hl636/rds/rds-mobile-Gwxt74iudAM/Public_genDiversity_databases/SPMP/dRep/dereplicated_genomes.txt`
	#bin=`sed -n ${counter}p < /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_UHGG_ASSEMBLY/5.4.DREP/dereplicated_genomes.txt`
	#bin=`sed -n ${counter}p < /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_UHGG_SPMP_ASSEMBLY/5.4.DREP/dereplicated_genomes.txt`
	#cat ${bin} | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/concat_metawrapIndAssembly_draft2.fa
	#cat ${bin} | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/concat_metawrapIndAssembly_draft2_sANI98.fa
	#cat ${bin} | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/concat_metawrapCoAssembly_draft1.fa
	#cat ${bin} | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/concat_metawrapCoAssembly_sANI98_draft2.fa
	#cat ${bin} | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/concat_metawrapHybridAssembly_Comp50Cont5_draft2.fa
	#cat ${bin} | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/concat_metawrapHybridAssembly_Comp70Cont5_draft1.fa
	#cat ${bin} | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/concat_metawrapHybridAssembly_sANI98_Comp50Cont5_draft2.fa
	#cat ${bin} | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/SPMPPolishedMAGs.fa
	cat ${bin} | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/SPMPPolishedMAGs_sANI95.fa
	#cat ${bin} | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/concat_metawrapHybrid_UHGG_Assembly_Comp50Cont5_draft2.fa
	#cat ${bin} | awk '/^>/{print ">MGYG""'$j'""_" ++i; next}{print}'  >> ${work_dir}/concat_metawrapHybrid_UHGG_SPMP_Assembly_Comp50Cont5_draft2.fa
	echo $counter
	counter=$((counter+1))
done

cd ${work_dir}
