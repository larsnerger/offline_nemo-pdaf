!$Id$
!BOP
!
! !ROUTINE: init_ens_offline --- Initialize ensemble for SEIK in offline mode
!
! !INTERFACE:
SUBROUTINE init_ens_offline(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize an ensemble of dim\_ens states.
! For the offline mode, the ensemble will be
! typically read-in from files.
!
! The routine is called by all filter processes and 
! initializes the ensemble for the PE-local domain.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
  use netcdf
  USE mod_assimilation_pdaf, &
       ONLY: program_mode, ensfile_type, timeDA, &
       flate, genEnsMeanYearly, nyears, GaussTransf, trafoConst, &
       flateZ, flateTOP, flateBOT, nLevFB, nLevFE
  use mod_io_pdaf, &
       only : read_state_mv, read_ens, gen_ens_mv, read_trajectory_mv, &
       write_state_ens, write_ens_files, write_ens_states, write_ens_fields, &
       gen_ensMeanYearly, &
       read_ens_dim_ens_files, read_ens_dims_dim_ens_files, &
       flate_depth, gen_ensFlateZ, nfiles, ntimec, &
       path_state, file_state_date1, file_state_date2, path_ens, file_ens, &
       ens_filelist, path_covar, file_covar, coupling_nemo
  USE mod_parallel_pdaf, &
       ONLY: mype=>mype_filter

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for SEIK. 
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_filter_init    (as U_ens_init)
!EOP

! *** local variables ***
  INTEGER :: i, member           ! Counters
  REAL :: invdim_ensm1                 ! Inverse of ensemble size minus 1
  CHARACTER(len=200):: titleVar
  REAL :: incrTime,startEnsTime,endEnsTime
  INTEGER, ALLOCATABLE :: hist_ens(:)
  REAL              :: delta,skewness,kurtosis
  INTEGER           :: status_hist,status_ensstats 
  INTEGER           :: ios,elemDiagn,k
  CHARACTER(len=20) :: rankfile
  REAL, ALLOCATABLE :: flate_z(:)

  integer(4) :: dim_state_ergom
  integer(4) :: dim_state_ensfile
  integer(4) :: dim_ens_ensfile


! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype==0) THEN
     WRITE (*, '(/a, 4x, a)') 'NEMO-PDAF', 'Initialize state ensemble'
     WRITE (*, '(a, 6x, a, i5)') 'NEMO-PDAF', '--- Ensemble size:  ', dim_ens
  END IF


! ***************************
! *** Check ensemble size ***
! ***************************
  
  ntimec = 1

!LN TODO: Move this to the actual ensemble reading
if (1==2) then
  if (ensfile_type == 3) then
     ! Determine ensemble size from number of files present
     call read_ens_dims_dim_ens_files(path_ens,dim_state_ensfile,dim_ens_ensfile)
  endif
endif


! ********************************
! *** Read ensemble from files ***
! ********************************

  ! Read ensemble central state vector state_p
  if (trim(program_mode)=='assim') then

     IF (mype==0) &
          write (*,'(a,4x,a)') 'NEMO-PDAF', '--- Read central model state of the ensemble'

     IF (coupling_nemo=='incr') THEN
        CALL read_state_mv(path_state, file_state_date1, file_state_date2, dim_p, timeDA, &
             coupling_nemo, state_p)
     ELSE
        CALL read_state_mv(path_state, file_state_date1, file_state_date2, dim_p, 1, &
             coupling_nemo, state_p)
     END IF
  endif

  IF (flateZ) THEN
     ALLOCATE (flate_z(dim_p))
     CALL flate_depth(dim_p,nLevFB,nLevFE,flateTOP,flateBOT,flate_z)
  ENDIF

  ! Read ensemble perturbations
  if (ensfile_type == 1) then

     IF (mype==0) &
          write (*,'(a,4x,a)') 'NEMO-PDAF','--- Initialize ensemble from single ensemble file'

     call read_ens(trim(path_ens)//trim(file_ens), dim_p, dim_ens, ens_p)

  else if  (ensfile_type == 2) then
     ! This option initializes ensemble perturbations by reading from a
     ! set of NEMO/ERGOM output files and subtracing the mean of these

     IF (.NOT. genEnsMeanYearly) THEN
        IF (.NOT. flateZ) THEN
           if (mype==0) write (*,'(a,4x,a)') &
                'NEMO-PDAF','--- Initialize ensemble from set of NEMO/ERGOM output files'

           CALL gen_ens_mv(flate, ens_filelist, path_state, dim_p, dim_ens, ens_p)
        ELSE
!            call gen_ensFlateZ(flate_z,varname,ens_filelist,path_ens,dim_p,dim_ens,GaussTransf,trafoConst,ens_p)
        ENDIF
     ELSE
!        call gen_ensMeanYearly(flate,varname,ens_filelist,path_ens,dim_p,dim_ens,nyears,GaussTransf,trafoConst,ens_p)!TO DO: domain decomp parallelisation (dim_p=dim_state at the moment)
     ENDIF
     IF (flateZ) DEALLOCATE (flate_z)
     IF (write_ens_states) call write_state_ens(path_ens, file_ens, dim_p, dim_ens, ens_p)

  else if (ensfile_type == 3) then

!     call read_ens_dim_ens_files(path_ens, dim_p, dim_ens, ens_p)

  else if (ensfile_type == 4) then

     ! This is the usual choice for program_mode='covar' 
     ! and could be used at the very beginning of a DA sequence

     if (mype==0) write (*,'(a, 4x,a)') 'NEMO-PDAF','--- Read model trajectory from NEMO/ERGOM output files'

     CALL read_trajectory_mv(flate, ens_filelist, path_state, dim_p, dim_ens, ens_p)

  else if (ensfile_type == 5) then

     ! Generate ensemble from covariance matrix by 2nd-order exact sampling
     ! This generates the full ensemble and sets the central state state_p

     if (mype==0) WRITE (*, '(a, 4x, a)') 'NEMO-PDAF','--- generate ensemble from covariance matrix'

     CALL gen_ens_from_cov(path_covar, file_covar, dim_p, dim_ens, state_p, ens_p)

  endif


  ! Add ensemble central state to perturbations
  if (trim(program_mode)=='assim' .and. ensfile_type /= 5) then
     write (*,'(a,4x,a)') 'NEMO-PDAF', '--- Add central model state to ensemble perturbations'

     DO member = 1, dim_ens
!$OMP PARALLEL DO  
        DO i=1, dim_p
           ens_p(i, member) = ens_p(i, member) + state_p(i)
        END DO
!$OMP END PARALLEL DO
     END DO
  end if

  ! write ensemble of files holding ensemble of fields
  ! Careful: This routine applied the backward state transformation
  IF (write_ens_fields) call write_ens_files(path_ens,'ensembleField',dim_p,dim_ens,ens_p)


! ****************
! *** clean up ***
! ****************

END SUBROUTINE init_ens_offline


subroutine gen_ens_from_cov(path, name, dim_p, dim_ens, state_p, ens_p)

  use pdaf_interfaces_module, only: PDAF_SampleEns
  use mod_io_pdaf, only: read_eof_cov_mv, path_covar, file_covar
  use mod_parallel_pdaf, &
       only: mype=>mype_filter, abort_parallel

  implicit none

! *** Arguments ***
  character(len=*), intent(in) :: path         !< File path
  character(len=*), intent(in) :: name         !< File name
  integer, intent(in) :: dim_p                 !< dimension of local state vector
  integer, intent(in) :: dim_ens               !< ensemble size
  real, intent(inout) :: state_p(dim_p)        !< state vector (central state of ensemble)
  real, intent(inout) :: ens_p(dim_p, dim_ens) !< ensemble vector

! *** Local variables ***
  integer :: rank
  integer :: status_pdaf
  real, allocatable :: eofV(:, :)
  real, allocatable :: svals(:)


! *****************************************
! *** Generate ensemble of model states ***
! *****************************************

  ! *** Rank of matrix is ensemble size minus one
  rank = dim_ens - 1

  ! allocate memory for temporary fields
  ALLOCATE(eofV(dim_p, rank))
  ALLOCATE(svals(rank))

  ! get eigenvalue and eigenvectors from file
  call read_eof_cov_mv(path_covar, file_covar, dim_p, rank, eofV, svals)

  ! *** Generate full ensemble on filter-PE 0 ***
  WRITE (*, '(a, 4x, a)') 'NEMO-PDAF','--- generate ensemble using PDAF_SampleEns'
  WRITE (*, '(a, 6x, a)') &
       'NEMO-PDAF','--- 2nd order exact sampling'
  WRITE (*, '(a, 6x, a, i5)') 'NEMO-PDAF','--- number of EOFs: ', rank

  ! Use PDAF routine to generate ensemble from covariance matrix
  CALL PDAF_SampleEns(dim_p, dim_ens, eofV, svals, state_p, ens_p, 1, status_pdaf)

  if (status_pdaf /= 0) then
     write (*, '(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in sample ensemble of PDAF - stopping! (PE ', mype, ')'
     call abort_parallel()
  end if

! ****************
! *** clean up ***
! ****************
  DEALLOCATE(svals, eofV)

end subroutine gen_ens_from_cov
