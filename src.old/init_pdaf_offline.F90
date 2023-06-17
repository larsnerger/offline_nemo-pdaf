!$Id$
!>  Interface routine to call initialization of PDAF
!!
!! This routine collects the initialization of variables for PDAF.
!! In addition, the initialization routine PDAF_init is called
!! such that the internal initialization of PDAF is performed.
!! This variant is for the offline mode of PDAF.
!!
!! __Revision history:__
!!
!!
subroutine init_pdaf()

  use netcdf
  use mod_parallel_pdaf, &     ! Parallelization variables
       only: mype=>mype_world, npes=>npes_world, n_modeltasks, task_id, &
       COMM_model, COMM_filter, COMM_couple, filterpe, abort_parallel
  use mod_assimilation_pdaf, & ! Variables for assimilation
       only: iday, dim_state_p, screen, filtertype, subtype, dim_ens, &
       incremental, covartype, type_forget, forget, &
       locweight, type_trans, type_sqrt, &
       ensfile_type, timeDA, &
       shiftObsInWet, flate, genEnsMeanYearly, nyears, GaussTransf, trafoConst, &
       EnsDiagnos, flateZ,flateTOP,flateBOT,nLevFB,nLevFE, &
       deg2rad, program_mode
  use mod_nemo_pdaf, &
       only: jpiglo, jpjglo, jpk, nwet, nwet3d, use_wet_state, &
       lon1, lat1, nldi, nlei, nldj, nlej
  use mod_io_pdaf, & 
       only: add_slash, write_var_time, write_ens_states, write_ens_fields, &
       saveState, coupling_nemo, saveIncr, do_deflate, &
       file_PDAF_state, file_PDAF_variance, file_PDAF_incr, ens_filelist, &
       path_state, file_state_date1, file_state_date2, path_ens, file_ens, &
       path_covar, file_covar, path_restart, file_restart, &
       startEnsTime, endEnsTime, incrTime
  use mod_statevector_pdaf, &
       only: setup_state
  use obs_prof_pdafomi, &
       only: assim_prof, rms_obs_prof, lradius_prof, sradius_prof, &
       path_prof, file_prof, variable_prof
  use obs_sst_cmems_pdafomi, &
      only: assim_sst_cmems, path_sst_cmems, file_sst_cmems, rms_obs_sst_cmems, &
            lradius_sst_cmems, sradius_sst_cmems, mode_sst_cmems, dist_sst_cmems
  use PDAFomi, &
       only: PDAFomi_set_domain_limits

       
  implicit none

! *** Local variables ***
  integer :: i                 ! Counter
  integer :: filter_param_i(7) ! Integer parameter array for filter
  real    :: filter_param_r(2) ! Real parameter array for filter
  integer :: status_pdaf       ! PDAF status flag
  integer :: dim2d, dim3d      ! Dimension of 2D and 3D field in state vector
  real :: lim_coords(2,2)      ! Limiting coordinates of sub-domain
  
  logical :: printconfig = .true.  ! Print information on all configuration parameters
    
  ! Namelist for PDAF settings
  namelist /pdaf_nml/ program_mode, &
       iday, filtertype, dim_ens, timeDA, screen, printconfig,   &
       forget, type_forget, locweight, &
       path_state,file_state_date1,  file_state_date2, &
       path_ens, file_ens, ensfile_type, path_covar, file_covar,  &
       file_PDAF_state, file_PDAF_variance, file_PDAF_incr,&
       ens_filelist, saveState, saveIncr, do_deflate, coupling_nemo, path_restart, file_restart,&
       shiftObsInWet, write_ens_fields, write_ens_states, flate, genEnsMeanYearly, nyears, GaussTransf, &
       EnsDiagnos,flateZ,flateTOP,flateBOT,nLevFE,nLEVFB, write_var_time, &
       assim_prof, rms_obs_prof, lradius_prof, sradius_prof, &
       path_prof, file_prof, variable_prof, &
       assim_sst_cmems, path_sst_cmems, file_sst_cmems, rms_obs_sst_cmems, lradius_sst_cmems,  &
       sradius_sst_cmems, mode_sst_cmems, dist_sst_cmems

  ! Namelist for ERGOM-PDAF settings
  namelist /ergompdaf_nml/ startEnsTime, endEnsTime, incrTime 
 
  ! External subroutines
  external :: init_ens_offline  ! Ensemble initialization
  

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  if (mype == 0) then
     write (*,'(/1x,a)') 'INITIALIZE PDAF - OFFLINE MODE'
  end if


! **********************************************************
! ***   CONTROL OF PDAF - used in call to PDAF_init      ***
! **********************************************************

  program_mode = 'assim'   ! Mode of the program: 
                   !   'assim' to perform analysis step 
                   !   'covar' for EOF decomposition to generate covariance matrix file

! *** IO options ***
  screen      = 2  ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
  filtertype = 7    ! Type of filter
                    !   (1) SEIK
                    !   (2) EnKF
                    !   (3) LSEIK
                    !   (4) ETKF
                    !   (5) LETKF
                    !   (6) ESTKF
                    !   (7) LESTKF
  dim_ens = 9       ! Size of ensemble for all ensemble filters
                    ! Number of EOFs to be used for SEEK
  subtype = 5       ! (5) Offline mode
  type_trans = 0    ! Type of ensemble transformation
                    !   SEIK/LSEIK and ESTKF/LESTKF:
                    !     (0) use deterministic omega
                    !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   ETKF/LETKF:
                    !     (0) use deterministic symmetric transformation
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
  type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                    !   (0) fixed
                    !   (1) global adaptive
                    !   (2) local adaptive for LSEIK/LETKF/LESTKF
  forget  = 1.0     ! Forgetting factor
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
  covartype = 1     ! Definition of factor in covar. matrix used in SEIK
                    !   (0) for dim_ens^-1 (old SEIK)
                    !   (1) for (dim_ens-1)^-1 (real ensemble covariance matrix)
                    !   This parameter has also to be set internally in PDAF_init.


! *********************************************************************
! ***   Settings for analysis steps  - used in call-back routines   ***
! *********************************************************************

! *** specifications for observations ***

  ! Which observations to assimilate
  assim_prof = .false.        ! Whether to assimilate profile data
  assim_sst_cmems = .false.   ! Whether to assimilate SST data from CMEMS

  ! Settings for profile observations
  rms_obs_prof = 0.8    ! Observation error stddev for profile data
  lradius_prof = 20.0  ! Radius in grid points for observation domain in local filters
  sradius_prof = lradius_prof  ! Support range for 5th-order polynomial
                    ! or range for 1/e for exponential weighting
  path_prof = '/cmems_archive/bm1302b/OBS/Batch/SHARK/d2019/'
  file_prof = '20191216shark_OXY.dat'

  ! Settings for CMEMS satellite SST
  rms_obs_sst_cmems = 0.8    ! Observation error stddev for SST data from CMEMS
  mode_sst_cmems = 1         ! Observation mode for SST_CMEMS: 
                             !  (0) linear interpolation to observation grid
                             !  (1) super-obbing: average 4 observation values
  lradius_sst_cmems = 10000.0_8 ! Radius in km for lon/lat (or in grid points)
  sradius_sst_cmems = lradius_sst_cmems

  
! *** Localization settings
  locweight = 2     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance

! *** File settings
  
  ensfile_type = 2  ! (1) for using ens file
                    ! (2) for using a set of single output files
                    ! (3) for using dim_ens ens files (ens(dim_p)) e.g. from SMHI toolbox
                    ! (4) for reading in trajectory for generating covariance matrix
                    ! (5) for init using covariance matrix file
  write_ens_states = .false. ! Whether to write ensemble file after generation for ensfile_type=2
  write_ens_fields = .false. ! Whether to write set of files holding ensemble fields
  timeDA = 1

  flate = 1.0 !Factor to inflate and deflate ensemble. inflation: flate > 1.0, deflation: flate < 1.0
  flateZ = .false. !if true, depth dependend flate factor is used.
  flateTOP = 1.0   !flate factor at top of domain
  flateBOT = 1.0   !flate factor at bottom of domain
  nLevFE = 30      !level from above at which flate value changes from linear interp. to flateBOT 
  nLevFB = 10      !level from above at which flate value changes from flateTOP to linear interp. between flate                   !TOP and flateBOT  

  genEnsMeanYearly = .false. !generates ensembles with mean of each year
  nyears= 4 !number of years taken for ensemble generation (only needed for genEnsMeanYearly=.TRUE.)

  GaussTransf = 0   ! (0) no transformation
                    ! (1) log with basis 10 (y=log_10 (x+A))
                    ! (2) ln (y=ln (x+A)) 
                    ! (3) box-cox

  trafoConst = 71.0D0 !TO DO: read it in 
  
  EnsDiagnos = .false. 


! *** Read namelist file for PDAF if the list is there ***
  open (500,file='pdaf.nml')
  read (500,NML=pdaf_nml)
  close (500)

  open (500,file='ergompdaf.nml')
  read (500,NML=ergompdaf_nml)
  close (500)


! *** Add trailing slash to paths ***
  call add_slash(path_prof)
  call add_slash(path_sst_cmems)
  call add_slash(path_state)
  call add_slash(path_ens)

! *** For generating the covariance matrix reset ensemble mode
  if (trim(program_mode)=='covar') ensfile_type = 4


! *** Activate netcdf deflation for single process runs
  if (npes==1) do_deflate = .true.

! *** Print configuration variables ***
  showconf: if (printconfig .and. mype==0) then

     write (*,'(/1x,a)') '-- Overview of PDAF configuration --'
     write (*,'(3x,a)') 'PDAF [pdaf_nml]:'
     write (*,'(5x,a,a)')     'program_mode       ', program_mode     
     write (*,'(5x,a,i10)')   'filtertype  ', filtertype
     write (*,'(5x,a,i10)')   'subtype     ', subtype
     write (*,'(5x,a,i10)')   'dim_ens     ', dim_ens
     write (*,'(5x,a,i10)')   'screen      ', screen
     write (*,'(5x,a,f10.2)') 'forget        ', forget
     write (*,'(5x,a,i10)')   'locweight   ', locweight
     write (*,'(5x,a,i10)')   'ensfile_type', ensfile_type
     write (*,'(5x,a,l)')     'assim_prof         ', assim_prof
     if (assim_prof) then
        write (*,'(6x,a,es10.2)')'lradius_prof      ', lradius_prof
        write (*,'(6x,a,es10.2)')'sradius_prof      ', lradius_prof
        write (*,'(6x,a,es10.2)')'rms_obs_prof      ', rms_obs_prof
        write (*,'(6x,a,a)')     'path_prof     ', trim(path_prof)
        write (*,'(6x,a,a)')     'file_prof     ', trim(file_prof)
     end if
     write (*,'(5x,a,l)')        'assim_sst_cmems    ', assim_sst_cmems
     if (assim_sst_cmems) then
        write (*,'(6x,a,es10.2)')'lradius_sst_cmems  ', lradius_sst_cmems
        write (*,'(6x,a,es10.2)')'sradius_sst_cmems  ', lradius_sst_cmems
        write (*,'(6x,a,es10.2)')'rms_obs_sst_cmems  ', rms_obs_sst_cmems
        write (*,'(6x,a,a)')     'path_sst_cmems    ', trim(path_sst_cmems)
        write (*,'(6x,a,a)')     'file_sst_cmems    ', trim(file_sst_cmems)
     end if
     write (*,'(1x,a/)') '-- End of PDAF configuration overview --'
     
  end if showconf  

  !!here trafo of rms_obs
  select case (GaussTransf)
  case(0)
     write(*,*) 'No Transformation of rms_obs_prof'
  case(1)
     write(*,*) 'use log basis 10 transformation of rms_obs_prof'
     rms_obs_prof = sqrt((log10(rms_obs_prof))**2.0D0)
  case(2)
     write(*,*) 'use ln transformation of rms_obs_prof'
     rms_obs_prof = sqrt((log(rms_obs_prof))**2.0D0)
  case(3)
     write(*,*) 'no transformation- box cox still needs to be implemented'
  case DEFAULT
     write(*,*) 'No Transformation of rms_obs_prof'
  end select


! ***********************************
! *** Some optional functionality ***
! ***********************************

! *** Parse command line options   ***
! *** This is optional, but useful ***

  call init_pdaf_parse()


! *** Initial Screen output ***
! *** This is optional      ***

  if (mype == 0) call init_pdaf_info()


! ************************************************
! *** Specify state vector and state dimension ***
! ************************************************

  call setup_state(dim_state_p)
  

! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! *****************************************************

  filter_param_i(1) = dim_state_p ! State dimension
  filter_param_i(2) = dim_ens     ! Size of ensemble
  filter_param_i(3) = 0           ! Smoother lag (not implemented here)
  filter_param_i(4) = incremental ! Whether to perform incremental analysis
  filter_param_i(5) = type_forget ! Type of forgetting factor
  filter_param_i(6) = type_trans  ! Type of ensemble transformation
  filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
  filter_param_r(1) = forget      ! Forgetting factor
     
  call PDAF_init(filtertype, subtype, 0, &
       filter_param_i, 7,&
       filter_param_r, 2, &
       COMM_model, COMM_filter, COMM_couple, &
       task_id, n_modeltasks, filterpe, init_ens_offline, &
       screen, status_pdaf)


! *** Check whether initialization of PDAF was successful ***
  if (status_pdaf /= 0) then
     write (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype,')'
     call abort_parallel()
  end if


! ******************************************************************
! *** Specify domain limits to limit observations to sub-domains ***
! ******************************************************************

  if (trim(dist_sst_cmems) == 'gp') then
     lim_coords(1,1) = nldi
     lim_coords(1,2) = nlei
     lim_coords(2,1) = nlej
     lim_coords(2,2) = nldj
  else
     lim_coords(1,1) = lon1(nldi)*deg2rad
     lim_coords(1,2) = lon1(nlei)*deg2rad
     lim_coords(2,1) = lat1(nlej)*deg2rad
     lim_coords(2,2) = lat1(nldj)*deg2rad
  end if

  call PDAFomi_set_domain_limits(lim_coords)

end subroutine init_pdaf
