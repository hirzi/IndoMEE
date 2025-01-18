#!/bin/bash

# Define working directories and variables
work_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/hybrid_assembly"
checkm2_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.2.CHECKM2"
dirs=(/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered)
prefixes=(IND CO)

### Collate checkm2 results
echo "######### Collating checkm2 results #########"
for i in 0 1; do
	prefix=${prefixes[$i]}
	cd ${checkm2_dir}/${prefix}
	echo ${i} ${prefix}
	if [ -f quality_report.tsv ]; then
   		if [[ ${i} == 0 ]]
		then
			cat quality_report.tsv > ${checkm2_dir}/checkm2_summary.txt
		else
			tail -n+2 quality_report.tsv >> ${checkm2_dir}/checkm2_summary.txt
		fi
	fi
done
# Remove unbinned contigs
grep -v "unbinned" ${checkm2_dir}/checkm2_summary.txt > ${checkm2_dir}/temp 
mv ${checkm2_dir}/temp ${checkm2_dir}/checkm2_summary.txt


### In some cases, the length of the checkm2 output may be shorter/less than the number of input bins.
# To check the length of the checkm2 output vs the number of input bins.
#for i in 0 1; do
#	dir=${dirs[$i]}
#	prefix=${prefixes[$i]}
#	cd ${dir}
#	bin_ori=$(ls | wc -l)
#	temp=$(cat ${checkm2_dir}/${i}/quality_report.tsv | wc -l)
#	bin_checkm=$(( $temp - 1 ))
#	echo ${i} ${bin_ori} ${bin_checkm}
#done
### Discrepancy between the length of the checkm2 output and the number of input bins can arise because checkm2 skips bins if no proteins could be inferred (e.g. WARNING: Skipping protein file x_bin.n.faa as it was empty.)
# To ensure our genomeInfo file for dRep matches the set of input bins, we'll find the failed bins and append them to the output
echo "######### Finding failed checkm2 runs #########"
if [ -f ${checkm2_dir}/checkm2_failed_bins ]; then
    rm ${checkm2_dir}/checkm2_failed_bins
fi
for i in 0 1; do
	dir=${dirs[$i]}
	prefix=${prefixes[$i]}
	cd ${dir}
	# Check that directory is not empty
	if [ -z "$(ls -A .)" ]; then
   		:
	else
   		ls -1 | sed -e 's/\.fa$//' | sort -V > ${checkm2_dir}/original_bins
		tail -n+2 ${checkm2_dir}/${prefix}/quality_report.tsv | cut -f1 | sort -V > ${checkm2_dir}/checkm2_successful_bins
		diff ${checkm2_dir}/original_bins ${checkm2_dir}/checkm2_successful_bins | grep "<" | cut -c 3- >> ${checkm2_dir}/checkm2_failed_bins
	fi
done
rm ${checkm2_dir}/original_bins ${checkm2_dir}/checkm2_successful_bins

### We pad the difference between the checkm2 completed runs and the total number of input bins by appending the failed bins.
# And format in accordance with dRep requirements. This will be the GENOMEINFO file to be called with the --genomeInfo GENOMEINFO flag.
# For reference, see: https://github.com/MrOlm/drep/issues/206 and https://drep.readthedocs.io/en/master/advanced_use.html) and with -d/--debug flag (i.e. debug flag, which allows drep to continue where it left off/last terminated). 
echo "######### Writing summary genomeInfo file for dRep #########"
cat ${checkm2_dir}/checkm2_summary.txt | awk 'NF{print $1 ".fa"}' > ${checkm2_dir}/first_col
cat ${checkm2_dir}/checkm2_summary.txt | cut -f2-3 > ${checkm2_dir}/second_col
cat ${checkm2_dir}/checkm2_failed_bins | awk 'NF{print $1 ".fa"}' | grep -v "unbinned" >> ${checkm2_dir}/first_col
n_failed_bins=$(cat ${checkm2_dir}/checkm2_failed_bins | grep -v "unbinned" | wc -l)
for i in $(seq 1 ${n_failed_bins}); do echo -e "0\t100"; done >> ${checkm2_dir}/second_col
paste ${checkm2_dir}/first_col ${checkm2_dir}/second_col > ${checkm2_dir}/temp
cat ${checkm2_dir}/temp | tr "\\t" "," > ${checkm2_dir}/temp2
echo "genome,completeness,contamination" > ${checkm2_dir}/checkm2_genInfo.csv
tail -n+2  ${checkm2_dir}/temp2 >> ${checkm2_dir}/checkm2_genInfo.csv
rm ${checkm2_dir}/temp ${checkm2_dir}/temp2 ${checkm2_dir}/first_col ${checkm2_dir}/second_col
