#!/bin/bash

#SBATCH -A JACOBS-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH --job-name=collateBins
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=./logs/collateBins_4dRep.%A.%a.log
#SBATCH --error=./errs/collateBins__4dRep.%A.%a.err

# Define working directories
working_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/hybrid_assembly"
output_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.4.DREP"
# Array of final assemblies (individual + co-assemblies)
dirs=(/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered)
prefixes=(IND CO)

# Remove files if already exist
rm ${output_dir}/allRefinedBins.txt

# Add sample ID prefix and collate bins
echo "#################################   Collating refined bins   #################################"
for i in 0 1; do
	dir=${dirs[$i]}
	prefix=${prefixes[$i]}
	echo ${dir} ${prefix}
	cd ${dir}
	for f in *.fa ; do
		if [[ ${f} != ${prefix}* ]]
		then
			mv -- "${f}" "${prefix}_${f}"
		fi
	done
	if [ -z "$(ls -A .)" ]; then
   		:
	else
   		ls -d $PWD/* >> ${output_dir}/allRefinedBins.txt
	fi
done
# Remove unbinned contigs
grep -v "unbinned.fa" ${output_dir}/allRefinedBins.txt > ${output_dir}/temp 
mv ${output_dir}/temp ${output_dir}/allRefinedBins.txt

echo "#################################   DONE   #################################"
cd ${working_dir}
