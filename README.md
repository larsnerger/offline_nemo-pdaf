# NEMO-PDAF offline

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
