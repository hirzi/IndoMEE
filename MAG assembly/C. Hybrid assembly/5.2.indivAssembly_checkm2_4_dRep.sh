#!/bin/bash

#SBATCH -A JACOBS-SL2-CPU
#SBATCH -p icelake-himem
#SBATCH --job-name=checkm2
#SBATCH--array=1-2%2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --output=./logs/checkm2.%A.%a.log
#SBATCH --error=./errs/checkm2.%A.%a.err

###############################################################
#### Notes ###
# For deciding between -c, --cpus-per-task vs -n, --ntasks; see: https://stackoverflow.com/questions/51139711/hpc-cluster-select-the-number-of-cpus-and-threads-in-slurm-sbatch and https://support.ceci-hpc.be/doc/_contents/SubmittingJobs/SlurmFAQ.html#Q05
# -c, --cpus-per-task vs -n, --ntasks. Note that generally: -n cores, --ntasks=cores should be used for MPI jobs/processes (distributed memory) and --ntasks=1 --cpus-per-task=cores for OpenMP jobs/process (shared memory). See: https://scicomp.ethz.ch/wiki/LSF_to_Slurm_quick_reference#Parallel_job
# The --ntasks value refers to the number of tasks to be launched by SLURM only. This usually equates to the number of MPI tasks launched. Reduce this from nodes*76 if demanded by memory requirements, or if OMP_NUM_THREADS>1.

###############################################################

# Install checkm2
#git clone --recursive https://github.com/chklovski/checkm2.git
#cd checkm2
#conda env create -n checkm2 -f checkm2.yml

# Install checkm2 diamond database to selected path
#conda activate checkm2
#checkm2 database --download --path /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/databases/CheckM2_Diamond

# Test checkm2
#checkm2 testrun

###############################################################
#### Define job/run parameters and variables ###

# Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#mem_gb=$(( $SLURM_MEM_PER_NODE / 1000 ))

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
#conda activate checkm2
conda activate /rds/user/hl636/hpc-work/conda_mamba/envs/checkm2

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
idx_0based=$(($IDX - 1))
output_dir="/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/HYBRID_ASSEMBLY/5.2.CHECKM2"
# Array of final assemblies (individual + co-assemblies)
dirs=(/home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered /home/hl636/rds/rds-mobile-Gwxt74iudAM/Metagenomics_genomeBasedReference/metaWRAP/OUTPUT/CO_ASSEMBLY/5.5.FILTER_GUNC/dereplicated_genomes_gunc_finalFiltered)
prefixes=(IND CO)
dir=${dirs[${idx_0based}]}
prefix=${prefixes[${idx_0based}]}

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

# Run checkm2. For reference, see: https://github.com/chklovski/CheckM2
# Parallelise across populations
echo "#################################   Running CheckM2 on refined bins   #################################"
echo -e "\nExecuting command:\n==================\n$CMD\n"
CMD="checkm2 predict --threads $SLURM_CPUS_PER_TASK -x .fa --input ${dir} --output-directory ${output_dir}/${prefix}"
eval $CMD

###############################################################
