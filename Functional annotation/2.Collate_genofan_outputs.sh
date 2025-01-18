#!/bin/bash

# Define variables
#input_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP/Comp50Cont5/dereplicated_genomes
input_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP_sfastANI_98/Comp50Cont5/dereplicated_genomes
#outdir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/func_annotation/genofan/Results/HybridAssembly_Comp50Cont5
outdir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/func_annotation/genofan/Results/HybridAssembly_sANI98_Comp50Cont5
#outdir_summary=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/func_annotation/genofan/Results/HybridAssembly_Comp50Cont5_Summary
outdir_summary=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/func_annotation/genofan/Results/HybridAssembly_sANI98_Comp50Cont5_Summary
#mkdir ${outdir} ${outdir_summary}
#mv ${input_dir}/*annotations ${outdir}/
cd ${outdir}

# Format of output files
#amrfinder_results.tsv # 1 header line
#kegg_modules.tsv # 1 header line
#kegg_orthologs.tsv # no header; only number; get orthology from https://www.genome.jp/tools/kofamkoala/
#cazy_results.tsv # no header; get CAZy acronym meanings from http://www.cazy.org/spip.php?article20
#eggnog.emapper.annotations # 4 suppl. header lines + 1 mand. header line + 3 suppl. footer lines
#antismash/summary.tsv # no header; genome_name, cluster
#gutsmash/summary.tsv # no header; genome_name, cluster_name, cluster, cluster_class

# Collate output files
if [[ ! -f ${outdir_summary}/amrfinder_results_summary.tsv && ! -f ${outdir_summary}/kegg_modules_summary.tsv && ! -f ${outdir_summary}/eggnog_emapper_summary.tsv && ! -f ${outdir_summary}/gutsmash_results_summary.tsv && ! -f ${outdir_summary}/antismash_results_summary.tsv && ! -f ${outdir_summary}/kegg_orthologs_summary.tsv && ! -f ${outdir_summary}/cazy_results_summary.tsv ]]
then
	j=1
	for i in $(ls -1d *annotations | sort -V); do
		echo ${i} ${j}
	   	cd ${i}
	   	if [[ ${j} == 1 ]]
		then
			cat amrfinder_results.tsv > ${outdir_summary}/amrfinder_results_summary.tsv
			cat kegg_modules.tsv > ${outdir_summary}/kegg_modules_summary.tsv
			tail -n+5 eggnog.emapper.annotations | head -n-3 > ${outdir_summary}/eggnog_emapper_summary.tsv
			echo -e 'genome\tcluster_name\tcluster\tcluster_class' > ${outdir_summary}/gutsmash_results_summary.tsv
			echo -e 'genome\tcluster' > ${outdir_summary}/antismash_results_summary.tsv
		else
			tail -n+2 amrfinder_results.tsv >> ${outdir_summary}/amrfinder_results_summary.tsv
			tail -n+2 kegg_modules.tsv >> ${outdir_summary}/kegg_modules_summary.tsv
			tail -n+6 eggnog.emapper.annotations | head -n-3 >> ${outdir_summary}/eggnog_emapper_summary.tsv
		fi
		cat kegg_orthologs.tsv >> ${outdir_summary}/kegg_orthologs_summary.tsv
		cat cazy_results.tsv >> ${outdir_summary}/cazy_results_summary.tsv
		cat gutsmash/summary.tsv >> ${outdir_summary}/gutsmash_results_summary.tsv
		cat antismash/summary.tsv >> ${outdir_summary}/antismash_results_summary.tsv
		j=$((j + 1))
		cd ..
	done
fi

