#!/bin/bash

#SBATCH -A JACOBS-SL3-CPU
#SBATCH -p sapphire
#SBATCH --job-name=metawrap_assembly
#SBATCH --array=1-89%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=180000
#SBATCH --time=12:00:00
#SBATCH --output=./logs/2.metawrap_assembly.%A.%a.log
#SBATCH --error=./errs/2.metawrap_assembly.%A.%a.err

###############################################################
#### Notes ###
# For deciding between -c, --cpus-per-task vs -n, --ntasks; see: https://stackoverflow.com/questions/51139711/hpc-cluster-select-the-number-of-cpus-and-threads-in-slurm-sbatch and https://support.ceci-hpc.be/doc/_contents/SubmittingJobs/SlurmFAQ.html#Q05
# -c, --cpus-per-task vs -n, --ntasks. Note that generally: -n cores, --ntasks=cores should be used for MPI jobs/processes (distributed memory) and --ntasks=1 --cpus-per-task=cores for OpenMP jobs/process (shared memory). See: https://scicomp.ethz.ch/wiki/LSF_to_Slurm_quick_reference#Parallel_job
# The --ntasks value refers to the number of tasks to be launched by SLURM only. This usually equates to the number of MPI tasks launched. Reduce this from nodes*76 if demanded by memory requirements, or if OMP_NUM_THREADS>1.

###############################################################
#### Define job/run parameters and variables ###

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
mem_gb=$(( $SLURM_MEM_PER_NODE / 1000 ))

# Define job index for job arrays
IDX=$SLURM_ARRAY_TASK_ID

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
mamba activate /home/hl636/rds/hpc-work/conda_mamba/envs/metawrap_env

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
cleaned_reads_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/1.5.CLEAN_READS"
assembly_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/2.ASSEMBLY"
fastq_list="Test3_fastqs_filtered.list"

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
echo -e "\nExecuting command:\n==================\n$CMD\n"

# Run command
#eval $CMD
# Step 2. Assemble metagenomes
# Choose between assembling each sample seperately, or co-assembling across all samples (e.g. in case percentage remapping of reads to assembly is low).
# Choose between megahit (--megahit; much faster) or metaspades (--metaspades)
# For co-assembly approach, we'll need to concatenate the reads from all the samples:
#cat ${cleaned_reads_dir}/*_1.fastq > ${cleaned_reads_dir}/ALL_READS_1.fastq
#cat ${cleaned_reads_dir}/*_2.fastq > ${cleaned_reads_dir}/ALL_READS_2.fastq
#metawrap assembly -1 ${cleaned_reads_dir}/ALL_READS_1.fastq -2 ${cleaned_reads_dir}/ALL_READS_2.fastq -m $SLURM_MEM_PER_NODE -t $SLURM_CPUS_PER_TASK --megahit -o ${assembly_dir}
# Parallelise across samples
sample=`sed -n ${IDX}p < ${workdir}/${fastq_list}`
# Remember, set allows you to define the elements of your list as variables, according to their order. E.g. for a two column file, the first column variable can be called with ${1} and the second with ${2}
set -- $sample
sample_prefix=$(basename "${1}" _1.fastq.gz)
sample_prefix_1=$(basename "${1}" .fastq.gz)
sample_prefix_2=$(basename "${2}" .fastq.gz)
# Run assembly
#metawrap assembly -1 ${cleaned_reads_dir}/${sample_prefix_1}.fastq -2 ${cleaned_reads_dir}/${sample_prefix_2}.fastq -m ${mem_gb} -t $SLURM_CPUS_PER_TASK --megahit -o ${assembly_dir}/${sample_prefix}
metawrap assembly -1 ${cleaned_reads_dir}/${sample_prefix_1}.fastq.gz -2 ${cleaned_reads_dir}/${sample_prefix_2}.fastq.gz -m ${mem_gb} -t $SLURM_CPUS_PER_TASK --megahit -o ${assembly_dir}/${sample_prefix}

###############################################################