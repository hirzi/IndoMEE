#!/bin/bash

#SBATCH -A JACOBS-SL3-CPU
#SBATCH -p sapphire
#SBATCH --job-name=collateBins
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=./logs/collateBins_4dRep.%A.%a.log
#SBATCH --error=./errs/collateBins__4dRep.%A.%a.err

# Define working directories
working_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/co_assembly"
input_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.BIN_REFINEMENT"
output_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.4.DREP"
pops="pop.list"
#pops=(ASMDAI, BA0, BAPDW, BRUSBG, MALHUL, MALRPN, MALSUL, PBSHF) # pop order

# Remove files if already exist
rm ${output_dir}/allRefinedBins.txt

# Add sample ID prefix and collate bins
echo "#################################   Collating refined bins   #################################"
for i in $(cat ${working_dir}/${pops}); do 
	echo ${i}
	cd ${input_dir}/${i}/metawrap_50_5_bins
	for f in *.fa ; do
		if [[ ${f} != ${i}* ]]
		then
			mv -- "${f}" "${i}_${f}"
		fi
	done
	if [ -z "$(ls -A .)" ]; then
   		:
	else
   		ls -d $PWD/* >> ${output_dir}/allRefinedBins.txt
	fi	
	cd ${input_dir}
done
# Remove unbinned contigs
grep -v "unbinned.fa" ${output_dir}/allRefinedBins.txt > ${output_dir}/temp 
mv ${output_dir}/temp ${output_dir}/allRefinedBins.txt

echo "#################################   DONE   #################################"
cd ${working_dir}/
