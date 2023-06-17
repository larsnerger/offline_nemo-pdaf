!$Id: prepoststep_ens_offline.F90 1589 2015-06-12 11:57:58Z lnerger $
!>  Used-defined Pre/Poststep routine for PDAF
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in all ensemble filters.
!! 
!! The routine is called for global filters (e.g. ESTKF)
!! before the analysis and after the ensemble transformation.
!! For local filters (e.g. LESTKF) the routine is called
!! before and after the loop over all local analysis
!! domains.
!!
!! The routine provides full access to the state 
!! estimate and the state ensemble to the user.
!! Thus, user-controlled pre- and poststep 
!! operations can be performed here. For example 
!! the forecast and the analysis states and ensemble
!! covariance matrix can be analyzed, e.g. by 
!! computing the estimated variances. 
!! For the offline mode, this routine is the place
!! in which the writing of the analysis ensemble
!! can be performed.
!!
!! If a user considers to perform adjustments to the 
!! estimates (e.g. for balances), this routine is 
!! the right place for it.
!!
!! Implementation for the 2D offline example
!! without model parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code based on offline_1D
!! * Later revisions - see repository log
!!
subroutine prepoststep_ens_offline(step, dim_p, dim_ens, dim_ens_p, dim_obs, &
     state_p, Uinv, ens_p, flag)

  use mpi
  use netcdf
  use mod_assimilation_pdaf, &
       only: dim_state_p, GaussTransf, trafoConst,&
       EnsDiagnos
  use mod_statevector_pdaf, &
       only: n_fields, id, sfields
  use mod_io_pdaf, &
       only : read_state_mv, write_field_mv, write_increment, update_restart_mv, &
       write_var_time, saveState, saveIncr, coupling_nemo, &
       file_PDAF_state, file_PDAF_incr, file_PDAF_variance, &
       file_state_date1, file_state_date2, path_state, path_restart, file_restart, &
       startEnsTime, endEnsTime, incrTime
  use mod_nemo_pdaf, &
       only: nlvls=>jpk, nlons=>jpiglo, nlats=>jpjglo, &
       depths=>gdept_1d, lons=>glamt, lats=>gphit, tmp_4d
  use mod_aux_pdaf, & 
       only: var_limits_mv
  use mod_parallel_pdaf, &
       only: mype=>mype_filter, comm_filter, MPIerr
  use mod_memcount_pdaf, &
       only: memcount

  implicit none

! *** Arguments ***
  integer, intent(in) :: step        !< Current time step (negative for call after forecast)
  integer, intent(in) :: dim_p       !< PE-local state dimension
  integer, intent(in) :: dim_ens     !< Size of state ensemble
  integer, intent(in) :: dim_ens_p   !< PE-local size of ensemble
  integer, intent(in) :: dim_obs     !< Dimension of observation vector
  real, intent(inout) :: state_p(dim_p) !< PE-local forecast/analysis state
  !< (The array 'state_p' is not generally not initialized in the case of SEIK.
  !< It can be used freely here.)
  real, intent(inout) :: Uinv(dim_ens-1, dim_ens-1) !< Inverse of matrix U
  real, intent(inout) :: ens_p(dim_p, dim_ens)      !< PE-local state ensemble
  integer, intent(in) :: flag        !< PDAF status flag


! *** local variables ***
  integer :: i, j, k, member, domain, cnt      ! counters
  integer, save :: allocflag = 0       ! Flag for memory counting
  logical, save :: firsttime = .true.  ! Routine is called for first time?
  real :: invdim_ens                   ! Inverse ensemble size
  real :: invdim_ensm1                 ! Inverse of ensemble size minus 1
  real, allocatable :: rmse_est_p(:)   ! PE-local estimated RMS errors (ensemble standard deviations)
  real, allocatable :: state_tmp(:)    ! temporary state vector to hold model state variances or increment
  integer,save :: writestep_var=1      ! Time index for file output of variance
  integer,save :: writestep_state=1    ! Time index for file output of state
  integer :: nsteps                    ! Number of steps written into file
  character(len=3) :: forana           ! String indicating forecast or analysis
!  real    :: limitOxyMin, limitOxyMax  ! Variable limits
  character(len=200) :: titleState, titleVar   ! Strings for file titles
  integer, allocatable :: dimfield_p(:) ! Local field dimensions
  integer, allocatable :: dimfield(:)  ! Global field dimensions
  real, allocatable :: rmse_est(:)     ! Global estimated RMS errors (ensmeble standard deviations)
  integer :: id_incrfield              ! Id of field in state vector that is written into increment file

  integer, allocatable :: hist_ens(:)
  real              :: delta, skewness, kurtosis
  integer           :: status_hist, status_ensstats 
  integer           :: ios, elemDiagn
  character(len=20) :: rankfile

  
! **********************
! *** INITIALIZATION ***
! **********************

  if (firsttime) then
     if (mype == 0) write (*, '(a, 5x, a)') 'NEMO-PDAF', 'Analyze forecasted state ensemble'
     forana = 'for'
  else
     if (mype == 0)write (*, '(a, 5x, a)') 'NEMO-PDAF', 'Analyze and write assimilated state ensemble'
     forana = 'ana'
  end if

  ! Allocate fields
  allocate(state_tmp(dim_p))
  if (firsttime) call memcount(3,'r',dim_p)

  ! Initialize numbers
  invdim_ens    = 1.0_8 / real(dim_ens,8)  
  invdim_ensm1  = 1.0_8 / real(dim_ens - 1,8)


! **************************************************************
! *** Perform prepoststep for ensemble filter.               ***
! *** The state and error information is completely in the   ***
! *** ensemble.                                              ***
! **************************************************************

  ! *** Compute mean state

  state_p = 0.0
  do member = 1, dim_ens
!$OMP PARALLEL DO  
     do i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i, member)
     end do
  end do
!$OMP PARALLEL DO  
  do i = 1, dim_p
    state_p(i) = invdim_ens * state_p(i)
  enddo

  
  if (EnsDiagnos) then
     write(*,*) "*** Ensemble diagnostics post***"
     allocate(hist_ens(dim_ens+1))
     if (firsttime) call memcount(3,'r',dim_ens+1)

    ! Read file with elem number for diagnostics
     write(*,*) "*** Reading elem list for diagnostics:"
     write(*,*) 'ElemForEnsStats.txt'

     open (unit=154,file='ElemForEnsStats.txt',iostat=ios)
     if (ios /= 0) write(*,*) 'Could not open file ElemForEnsStats.txt'

     open (157, file='PostEnsDiagnostics.dat',iostat=ios)
       write(157,*)'elem  skewness  kurtosis'

     j=1
     hist_ens=0.0D0
     do
       read (154,*,iostat=ios) elemDiagn
       if (ios/=0) exit   
       call PDAF_diag_ensstats(dim_p, dim_ens,elemDiagn,state_p,ens_p,skewness,kurtosis,status_ensstats)
       write(157,*)elemDiagn,skewness,kurtosis
       call PDAF_diag_histogram(j,dim_p,dim_ens,elemDiagn,state_p,ens_p,hist_ens,delta,status_hist)
       j=j+1
     enddo
     close(154)
     rankfile=trim('PostRankHis_Diag.dat')
     print*,'rankfile',rankfile
     open (159, file=rankfile, iostat=ios)
         do k = 1, (dim_ens+1)
            write(159,'(I2,I10)')k,hist_ens(k)
         enddo
     close(159)


    ! ALLOCATE(hist_ens(dim_ens+1))
     hist_ens = 0.0D0
     !CALL PDAF_diag_histogram(1,dim_p,dim_ens,0,state_p,ens_p,hist_ens,delta,status_hist)
     !CALL PDAF_diag_ensstats(dim_p, dim_ens,0,state_p,ens_p,skewness,kurtosis,status_ensstats)
     call diag_histogram(1,dim_p,dim_ens,0,state_p,ens_p,hist_ens,delta,status_hist)
     call diag_ensstats(dim_p, dim_ens,0,state_p,ens_p,skewness,kurtosis,status_ensstats)
     !Write .dat file with j hist(j) delta skewness kurtosis
     open (158, file='PostRankHistogr.dat', iostat=ios)
         do j = 1, (dim_ens+1)
            write(158,'(I2,I10)')j,hist_ens(j)
         enddo
     close(158)
    ! OPEN (152, file='EnsDiagnostics.dat',status='new')
       !  WRITE(152,*)'skewness  kurtosis'
         write(157,*)'integral',skewness,kurtosis
     close(157)
     deallocate(hist_ens)
  endif


! LN: Limits are now set in setup_state and applied in var_limits_mv (in mod_aux)
  !here trafo of limit
!   select case (GaussTransf)
!      case(0)
!         write(*,*) 'No Transformation of bio limit'
!         if (varname == 'OXY') then
!            limitOxyMin = -450.0D0
!            limitOxyMax = 450.0D0
!         elseif (varname == 'NO3') then
!            limitOxyMin = 0.0D0
!            limitOxyMax = 8000.0D0 
!         elseif (varname == 'NH4') then
!            limitOxyMin = 0.0D0
!            limitOxyMax = 8000.0D0
!         elseif (varname == 'PO4') then
!            limitOxyMin = 0.0D0
!            limitOxyMax = 8000.0D0
!         else
!            limitOxyMin = -70.0D0
!            limitOxyMax = 450.0D0
!         endif
!      case(1)
!         write(*,*) 'use log basis 10 transformation of bio variable'
!         write(*,*) 'Attention: only implemented for oxygen variable'
!         limitOxyMin = log10(-70.0D0+trafoConst)
!         limitOxyMax = log10(450.0D0+trafoConst)
!      case(2)
!         write(*,*) 'use ln transformation of bio variable'
!         write(*,*) 'Attention: only implemented for oxygen variable'
!         limitOxyMin = log(-70.0D0+trafoConst)
!         limitOxyMax = log(450.0D0+trafoConst)
!      case(3)
!         write(*,*) 'no transformation- box cox still needs to be implemented'
!         limitOxyMin = -70.0D0
!         limitOxyMax = 450.0D0
!      case DEFAULT
!         write(*,*) 'No Transformation of bio variable'
!         limitOxyMin = -70.0D0
!         limitOxyMax = 450.0D0
!   end select

  ! Limits to each variables are applied according to the definitions in sfields
  ! The limits are also applied when applying the variable transformation for file output
  call var_limits_mv(state_p)


  ! *** Compute sampled variances ***
  state_tmp(:) = 0.0
  
  do member = 1, dim_ens
!$OMP PARALLEL DO  
     do j = 1, dim_p
        state_tmp(j) = state_tmp(j) &
             + (ens_p(j, member) - state_p(j)) &
             * (ens_p(j, member) - state_p(j))
     end do
  end do
!$OMP PARALLEL DO  
  do j = 1, dim_p
     state_tmp(j) = invdim_ensm1 * state_tmp(j)
  enddo
  

! ************************************************************
! *** Compute RMS errors according to sampled covar matrix ***
! ************************************************************

  allocate(rmse_est_p(n_fields))
  allocate(rmse_est(n_fields))
  allocate(dimfield_p(n_fields))
  allocate(dimfield(n_fields))
  rmse_est_p  = 0.0_8

  dimfield_p(:) = sfields(:)%dim
  CALL MPI_Allreduce(dimfield_p, dimfield, n_fields, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)

  ! total estimated mean-square error per field per process
  do j = 1, n_fields

     do i = 1+sfields(j)%off, sfields(j)%dim+sfields(j)%off
        rmse_est_p(j) = rmse_est_p(j) + state_tmp(i)
     enddo
     rmse_est_p(j) = rmse_est_p(j) / real(dimfield(j), 8)

  enddo

  ! Global sum of mean squared errors
  CALL MPI_Allreduce (rmse_est_p, rmse_est, n_fields, MPI_DOUBLE_PRECISION, MPI_SUM, &
       COMM_filter, MPIerr)

  ! Get global RMSE
  rmse_est = sqrt(rmse_est)


! *****************
! *** Screen IO ***
! *****************

  ! Output RMS errors given by sampled covar matrix
   if (mype == 0) then
      write (*, '(a,6x,a)') 'NEMO-PDAF', 'RMS errors according to sample variance'
      do i = 1, n_fields
         write (*,'(a,4x,a8,4x,a10,2x,es12.4)') &
              'NEMO-PDAF', 'RMSE-'//forana, trim(sfields(i)%variable), rmse_est(i)
      end do
   end if


! *******************
! *** File output ***
! *******************

   ! Write variance into nc file  
   if ((trim(write_var_time)=='fcst' .and. firsttime) .or. trim(write_var_time)=='both') then
      if (trim(write_var_time)=='fcst') then

         titleVar='Forecast ensemble variance'
         if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write variance before analysis step'
      else
         titleVar='Ensemble variance'
         if (mype == 0) write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write variance file'
      end if

      if (write_var_time=='both') then
         nsteps = 2
      else
         nsteps = 1
      end if

      call write_field_mv(state_tmp, dim_p, file_PDAF_variance, titleVar, &
           incrTime, nsteps, writestep_var)

      writestep_var = writestep_var + 1

   elseif (trim(write_var_time)=='ana' .and. (.not.firsttime)) then

      write (*,'(a,5x,a)') 'NEMO-PDAF', '--- Write variance after analysis step'
      titleVar='Analysis ensemble variance'

      call write_field_mv(state_tmp, dim_p, file_PDAF_variance, titleVar, incrTime, 2, 1)

   end if

   ! Write state into nc file
   if (saveState) then

      ! Store state in state_tmp to avoid changing state_p
      state_tmp = state_p

      ! Write state file for viewing
      if (forana=='for') then
         titleState = 'Forecast State'
         writestep_state = 1
      else
         titleState = 'Analysis State'
         writestep_state = 2
      end if

      ! Write separate files for forecast and analysis
!      call write_field_mv(state_tmp, dim_p, trim(file_PDAF_state)//'_'//forana//'.nc', titleState, incrTime, 1, 1)
      ! Write forecast and analysis into the same file
      call write_field_mv(state_tmp, dim_p, trim(file_PDAF_state)//'.nc', titleState, incrTime, 2, writestep_state)
   endif


   ! *** Operations after analysis ***

   notfirst: if (.not. firsttime) then

      if (mype==0) write (*,'(a,2x,a)') &
           'NEMO-PDAF', 'Read forecast model state to compute increment'

      ! Calculate assimilation increment

      call read_state_mv(path_state, file_state_date1, file_state_date2, dim_p, 1, &
           coupling_nemo, state_tmp)

      state_tmp=state_p-state_tmp


      ! Write increment to netCDF for input in ERGOM
      if (coupling_nemo=='incr' .or. saveIncr) then
    
         id_incrfield = id%temp

         call write_increment(state_tmp, dim_p, file_PDAF_incr, id_incrfield)

      endif

      ! Update restart file when coupling to NEMO is done through this file
      if(coupling_nemo=='rest') then

         call update_restart_mv(path_restart, file_restart, state_p, state_tmp)

      endif

   end if notfirst



! ********************
! *** finishing up ***
! ********************

   deallocate(state_tmp) 
   deallocate(rmse_est_p, rmse_est, dimfield_p, dimfield)

   firsttime = .false.

end subroutine prepoststep_ens_offline
