#!/bin/bash
#SBATCH --job-name=PDAFOFF
###SBATCH -p smp 
#SBATCH -p mpp
##SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=0:30:00
##SBATCH --mail-user=lars.nerger@awi.de
##SBATCH --mail-type=END
#
#

ulimit -s unlimited
export LD_LIBRARY_PATH=/global/AWIsoft/hdf5-parallel/1.12.1_intel2021/lib:${LD_LIBRARY_PATH}":"
echo $LD_LIBRARY_PATH

export OMP_NUM_THREADS=3
echo $OMP_NUM_THREADS

srun ./PDAF_offline 


