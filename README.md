# NEMO-PDAF offline

This code is aimed at a particular test case for the North and Baltic Seas performing a single analysis assimilating sea surface observations. WE provide further information on the particular test case further below.

## Preparation for compilation

Next to this NEMO-PDAF offline code you should get PDAF from https://github.com/PDAF/PDAF. 
You should put the PDAF directory with the name `PDAF` into the main directory of the offline code.

The compilation uses the usual build scheme of PDAF (see. e.g. https://pdaf.awi.de/trac/wiki/FirstSteps). 
Thus, you should choose a file with build definitions from `PDAF/make.arch/`. If there is no suitable file, 
you can create your own based on the existing files. (Note for JUWELS: For JUWELS, there is the file `setup/JUWELS_mpifort_impi.h` in the offline_nemo-pdaf code directory. This should be copied into `PDAF/make.arch/`. Then one can compile using the specification `PDAF_ARCH=JUWELS_mpifort_impi`)

You should compile with support for MPI and OpenMP. The compilation needs the BLAS and LAPACK libraries 
(or some scientific library which provides the functions of these libraries). Further, a NetCDF library is required.
In case of using the MPI parallelization, you need a parallel NetCDF library. 

Note for JUWELS: Here we use the configuration (as of July 2023)<br>
`module load Intel/2022.1.0`<br>
`module load ParaStationMPI/5.8.1-1`<br>
`module load  imkl/2022.1.0`<br>
`module load  netCDF-Fortran/4.6.0`

## Compilation

cd in to src/ and run `make`. This will compile both the PDAF library and the NEMO-PDAF offline program `PDAF_offline`.

## Running

### Preparare the run

You can run the program in `run/`. 

To preparare the run you should
* Place the observation file in `obs/`. The default name is `sst_multi_201801.nc` (multi-satellite SST from CMEMS).
* Place the input files (snapshots) from the free run into `out_free/`.
* Link the executable `PDAF_offline` from `src/` into `run/`
* Copy the files from `setup/` into `run/`.
* The run script `run_pdaf_1task.sl` is for runnig with SLURM (configured for JUWELS). This file is configured to use a single compute node with 48 processor cores to run with a single MPI task and 48 OpenMP threads. The files `run_pdaf_2tasks.sl` and `run_pdaf_4tasks.sl` are configured for 2 MPI tasks using 24 OpenMP threads each or 4 MPI tasks using 12 OpenMP threads each, respectively. The run script `run_pdaf_2node.sl` uses 2 MPI tasks with 48 OpenMP threads each. 

Now the job is ready to run. For running with SLURM, you can execute `sbatch run_pdaf_1task.sl`. The run produces the output files `PDAFstate.nc`, `PDAFVariance.nc` and `PDAFIncrement.nc`. On JUWELS with 48 OpenMP threads the job should take about 3 minutes to run. The job output shows detailed timing and memory information at the end. 

The job needs about 7 GB of RAM. The outputs have a size of about 6.6 GB.

## The test case

The test case is a realistic setup of moderate size (a similar case is currently run operational to produce forecasts for the Baltic Sea within the Copernicus Marine Service). The configuration is for using the shared-memory parallelization with OpenMP, but also moderate MPI-parallelization can be used (generally, one can also run with a high number of MPI tasks, but the code has no functionality to automatically generate the decomposition, so that one would need to specify sub-domains manually in a configuration file. We run this case in a variant in which PDAF is directly coupled into the ocean model, but this online-coupled case looks to complicated as a test case for porting to GPUs.) 

The test case performs a single analysis (assimilation) step that assimilates a field of sea surface temperature. Used is the ensemble Kalman filter LESTKF with an ensemble of 30 members. The ensemble is read before the analysis in `PDAF_init` and the ensemble mean and the ensemble variance before and after the assimialtion update are writting into netcdf files. Further an increment file is written. This could then be used with the NEMO ocean model to apply the increment (Here we omit the use of NEMO and focus on the analysis step of the data assimilation).

The motivation for porting to GPUs comes from the fact that more-and-more Earth system component models are ported to GPUs and these aim for exascale. For such models, also the data assimilation component needs GPU support. 

An important aspect of PDAF is that there is a generic part that is computed inside the library itself (mainly doing linear algebra) and case-specific routines, which are called by PDAF as call-back functions. The case-specific routines are those in `offline_nemo-pdaf/src`, while the generic are in the PDAF library `PDAF/src`. Porting to GPU should focus on the operations inside the calls to PDAF_init and PDAFomi_assimilate_local and here in particular the operation inside the PDAF library. Since the call-back routines are case-specific, they are usually coded by users based on templates provided with the PDAF code (in PDAF/templates/). Porting the routines of this test case would be useful as it provides hints on the porting also for other cases. 

Note, that the major execution time in the actual assimilation step inside `PDAFomi_assimilate_local` (in assimilation_pdaf_offline.F90) shifts from the call-back routine `init_dim_obs_pdafomi` (in callback_obs_pdafomi.F90) for a single MPI task to operations inside the PDAF library for a larger number of MPI tasks. The execution time for the initialization in PDAF_init, and here in particular the time tine `init_ens_pdaf` should be of less concern. (In the case that PDAF is directly coupling into a model, this would only be executed once at the very beginning of a data assimilation sequence). The major candidate for porting should be the calculations inside the analysis step (listed under `LESTKF analysis` in the timing output).


