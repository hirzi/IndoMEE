#!/bin/bash

#SBATCH -A JACOBS-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH --job-name=dRep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=./logs/dRep.%A.%a.log
#SBATCH --error=./errs/dRep.%A.%a.err

###############################################################
#### Notes ###
# For deciding between -c, --cpus-per-task vs -n, --ntasks; see: https://stackoverflow.com/questions/51139711/hpc-cluster-select-the-number-of-cpus-and-threads-in-slurm-sbatch and https://support.ceci-hpc.be/doc/_contents/SubmittingJobs/SlurmFAQ.html#Q05
# -c, --cpus-per-task vs -n, --ntasks. Note that generally: -n cores, --ntasks=cores should be used for MPI jobs/processes (distributed memory) and --ntasks=1 --cpus-per-task=cores for OpenMP jobs/process (shared memory). See: https://scicomp.ethz.ch/wiki/LSF_to_Slurm_quick_reference#Parallel_job
# The --ntasks value refers to the number of tasks to be launched by SLURM only. This usually equates to the number of MPI tasks launched. Reduce this from nodes*76 if demanded by memory requirements, or if OMP_NUM_THREADS>1.
# Note: blobology is very resource intensive, so either allocate a high number of cores, allocate a longer runtime (SL2), or subsample the run.
# --array=1-3%3

###############################################################
#### Define job/run parameters and variables ###

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#mem_gb=$(( $SLURM_MEM_PER_NODE / 1000 ))

# Define job index for job arrays
#IDX=$SLURM_ARRAY_TASK_ID

# Load modules
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment
module load miniconda/3

# Load environments.
# Note that conda/mamba init + source ~/.bashrc modifies your bashrc script, sometimes resulting in artifacts in the bashrc script (when running in job array). So 
eval "$(conda shell.bash hook)"
#source /rds/user/hl636/hpc-work/conda_mamba/etc/profile.d/conda.sh
#conda init
#source ~/.bashrc
conda activate /rds/user/hl636/hpc-work/conda_mamba
#mamba init
#source ~/.bashrc
source /rds/user/hl636/hpc-work/conda_mamba/etc/profile.d/mamba.sh
#mamba activate metawrap_env
#mamba activate metawrap_env_py3
#mamba activate /home/hl636/rds/hpc-work/conda_mamba/envs/metawrap_env
mamba activate /home/hl636/rds/hpc-work/conda_mamba/envs/metawrap_env_py3

#! Full path to application executable:
#application=

#! Run options for the application:
#options=

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"

#! Are you using OpenMP (NB this is unrelated to OpenMPI)? If so increase this safe value to no more than 76:
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$[${numnodes}*${mpi_tasks_per_node}]

#! The following variables define a sensible pinning strategy for Intel MPI tasks -
#! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets
#! Notes:
#! 1. These variables influence Intel MPI only.
#! 2. Domains are non-overlapping sets of cores which map 1-1 to MPI tasks.
#! 3. I_MPI_PIN_PROCESSOR_LIST is ignored if I_MPI_PIN_DOMAIN is set.
#! 4. If MPI tasks perform better when sharing caches/sockets, try I_MPI_PIN_ORDER=compact.

#! Uncomment one choice for CMD below (add mpirun/mpiexec options if necessary):

#! Choose this for a MPI code (possibly using OpenMP) using Intel MPI.
#CMD="mpirun -ppn $mpi_tasks_per_node -np $np $application $options"

#! Choose this for a pure shared-memory OpenMP parallel program on a single node:
#! (OMP_NUM_THREADS threads will be created):
#CMD="$application $options"

#! Choose this for a MPI code (possibly using OpenMP) using OpenMPI:
#CMD="mpirun -npernode $mpi_tasks_per_node -np $np $application $options"

###############################################################
#### RUN ###

# Change working directory
cd $workdir
echo -e "Changed directory to `pwd`.\n"

# Set other directories
working_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/co_assembly"
checkm2_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.2.CHECKM2"
input_bins=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.4.DREP/allRefinedBins.txt
genome_info=${checkm2_dir}/checkm2_genInfo.csv
output_dir=/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.4.DREP
#mkdir ${output_dir}

# Print job info
JOBID=$SLURM_JOB_ID
echo -e "JobID: $JOBID\n======"
echo "Time: `date`"
echo "Running on master node: `hostname`"
echo "Current directory: `pwd`"
if [ "$SLURM_JOB_NODELIST" ]; then
        echo -e "\nNodes allocated:\n================"
fi
echo -e "\nnumtasks=$numtasks, numnodes=$numnodes, mpi_tasks_per_node=$mpi_tasks_per_node (OMP_NUM_THREADS=$OMP_NUM_THREADS)"

# Run command
# Dereplicate bins (for individual assembly strategy). For reference, see: https://drep.readthedocs.io/en/latest/module_descriptions.html ; https://drep.readthedocs.io/en/latest/choosing_parameters.html
# -g # Can also input a text file  with paths to genomes
# When continuing run (via -debug; https://github.com/MrOlm/drep/issues/83), you may need to remove -g XXX
# When running for many genomes (i.e. > --primary_chunksize; with args --multiround_primary_clustering --primary_chunksize 5000), you may encounter an error (drep may freeze) when partitioning genomes into multiround groups. Error below:
# > Running pair-wise MASH clustering
# > Will split genomes into 2 groups for primary clustering (freezes here)
# The solution is either to increase primary chunksize to avoid multiround behaviour, or for large sets of genomes, you'll need to debug multiround_primary_clustering properly
#CMD="dRep dereplicate ${output_dir} -g ${input_bins} -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI --genomeInfo ${genome_info} -pa 0.9 -sa 0.95 -nc 0.1 -l 50000 -comp 75 -con 25 --multiround_primary_clustering --primary_chunksize 5000 --greedy_secondary_clustering --run_tertiary_clustering --debug" # default
#CMD="dRep dereplicate ${output_dir} -g ${input_bins} -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI --genomeInfo ${genome_info} -pa 0.9 -sa 0.95 -nc 0.1 -l 50000 -comp 70 -con 5 --multiround_primary_clustering --primary_chunksize 5000 --greedy_secondary_clustering --run_tertiary_clustering --debug" # match metawrap, with greedy seconday clustering (a little bit faster)
#CMD="dRep dereplicate ${output_dir} -g ${input_bins} -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI --genomeInfo ${genome_info} -pa 0.9 -sa 0.95 -nc 0.1 -l 50000 -comp 70 -con 5 --multiround_primary_clustering --primary_chunksize 5000 --run_tertiary_clustering --debug" # match metawrap
#CMD="dRep dereplicate ${output_dir} -g ${input_bins} -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI --genomeInfo ${genome_info} -pa 0.9 -sa 0.95 -nc 0.3 -l 50000 -comp 50 -con 5 --multiround_primary_clustering --primary_chunksize 5000 --run_tertiary_clustering --debug" # revised parameters
#dRep dereplicate ${output_dir} -g ${input_bins} -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI --genomeInfo ${genome_info} -pa 0.9 -sa 0.95 -nc 0.3 -l 50000 -comp 50 -con 5 --multiround_primary_clustering --primary_chunksize 10000 --run_tertiary_clustering --debug # revised parameters
dRep dereplicate ${output_dir} -g ${input_bins} -p $SLURM_CPUS_PER_TASK --S_algorithm fastANI --genomeInfo ${genome_info} -pa 0.9 -sa 0.98 -nc 0.3 -l 50000 -comp 50 -con 5 --multiround_primary_clustering --primary_chunksize 10000 --run_tertiary_clustering --debug # revised parameters
#eval $CMD
#echo -e "\nExecuting command:\n==================\n$CMD\n"

###############################################################