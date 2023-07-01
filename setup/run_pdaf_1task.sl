#!/bin/bash
#SBATCH --account=esmtst
##SBATCH --job-name=PDAFOFF
##SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --time=0:30:00
#SBATCH --mail-user=lars.nerger@awi.de
#SBATCH --mail-type=END
#SBATCH --partition=batch   
#
#
#source ${MODULESHOME}/init/bash
source /etc/profile.d/modules.sh

#module use /global/AWImodules/
#module purge
#module load intel.compiler/2021.3.0
#module load intel.mpi/2021.3.0
#module load netcdf-parallel/4.8.1_intel2021

ulimit -s unlimited
export LD_LIBRARY_PATH=/global/AWIsoft/hdf5-parallel/1.12.1_intel2021/lib:${LD_LIBRARY_PATH}":"
echo $LD_LIBRARY_PATH

wdir=`pwd`
echo ' '
echo 'Run directory: ' $wdir
export wdir

export OMP_PROC_BIND=TRUE
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
echo 'OMP_NUM_THREADS ' $OMP_NUM_THREADS

srun ./PDAF_offline 


