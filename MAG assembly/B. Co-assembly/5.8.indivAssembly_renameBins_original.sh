#!/bin/bash

#SBATCH -A JACOBS-SL3-CPU
#SBATCH -p sapphire
#SBATCH --job-name=renameBins_original
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=./logs/renameBins_original.%A.%a.log
#SBATCH --error=./errs/renameBins_original.%A.%a.err

# Define working directories
working_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/co_assembly"
input_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.BIN_REFINEMENT"
pops="pop.list"

# Remove sample ID prefix from bin names
echo "#################################   Renaming bins to original names  #################################"
for i in $(cat ${working_dir}/${pops}); do 
	echo ${i}
	cd ${input_dir}/${i}/metawrap_50_5_bins
	for filename in *.fa; do
		if [[ ${filename} == ${i}* ]]
		then	
			[ -f "$filename" ] || continue
			mv "$filename" "${filename//${i}_/}"
		fi
	done
	cd ${input_dir}
done

echo "#################################   DONE   #################################"
cd ${working_dir}/
