# NEMO-PDAF offline

This code is aimed at a particular test case for the North and Baltic Seas performing a single analysis assimilating sea surface observations.

## Preparation for compilation

Next to this NEMO-PDAF offline code you should get PDAF from https://github.com/PDAF/PDAF. 
You should put the PDAF directory with the name `PDAF` into the main directory of the offline code.

The compilation uses the usual build scheme of PDAF (see. e.g. https://pdaf.awi.de/trac/wiki/FirstSteps). 
Thus, you should choose a file with build definitions from `PDAF/make.arch/`. If there is no suitable file, 
you can create your own based on the existing files.

You should compile with support for MPI and OpenMP. The compilation needs the BLAS and LAPACK libraries 
(or some scientific library which provides the functions of these libraries). Further, a NetCDF library is required.
In case of using the MPI parallelization, you need a parallel NetCDF library. 

## Compilation

cd in to src/ and run `make`.

## Running

### Preparare the run

You can run the program in `run/`. 

To preparare the run you should
* Place the observation file in `obs/`. The default name is `sst_multi_201801.nc` (multi-satellite SST from CMEMS).
* Place the input files (snapshots) from the free run into `out_free/`.
* Link the executable `PDAF_offline` from `src/` into `run/`
* Copy the files from `setup/` into `run/`.
* The run script `run_pdaf_1pe.sl` is for runnig with SLURM. This file is configured to for a compute node with 128 processor cores to run with a single MPI task and 128 OpenMP threads. The file `run_pdaf_4pe.sl` is configures for 4 MPI tasks each using 32 OpenMP threads.

Now the job is ready to run. For running with SLURM, you can execute `sbatch run_pdaf_1pe.sl`. The run produces the output files `PDAFstate.nc`, `PDAFVariance.nc` and `PDAFIncrement.nc`. With 128 OpenMP threads the job should take about 2 minutes to run.

The job needs about 7 GB of RAM. The output have a size of about 750 MB.




