!$Id$
!BOP
!
! !ROUTINE: eofcovar - Control routine for generating covariance matrix
!
! !INTERFACE:
subroutine eofcovar()

! !DESCRIPTION:
! This routine performs a single analysis step for
! PDAF in offline mode using PDAF with domain-decomposition.
!
! The analysis is performed by calling a filter-specific 
! routine PDAF\_put\_state\_X.
!
! In this routine, the real names of most of the 
! user-supplied routines for PDAF are specified (see below).
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code by restructuring
! Later revisions - see svn log
!
! !USES:
  use PDAF_interfaces_module, &
       only: PDAF_set_ens_pointer
  use mod_parallel_pdaf, &    ! Parallelization
       only: mype_world, abort_parallel
  use mod_assimilation_pdaf, & ! Variables for assimilation
       only: do_mv, remove_mean, limit, hwindow, &
       dim_ens, dim_state_p
  use mod_statevector_pdaf, &
       only: n_fields, sfields
  use mod_io_pdaf, &
       only: write_covar_mv, path_covar, file_covar
  use mod_memcount_pdaf, &
       only: memcount
  
  implicit none

! local variables
  integer :: i, istep
  integer :: status    ! Status flag for filter routines
  integer :: maxtimes  
  real, pointer :: states(:,:)            ! Pointer to state trajectory (PDAF's ensemble array)
  real, allocatable :: meanstate(:)       ! time mean state vector
  real, allocatable :: tmp_ens(:,:) ! temporary: running mean states / left singular vectors
  integer :: rank
  real, allocatable :: svals(:) ! field of singular values
  real, allocatable :: stddev(:) ! stddev of field

  namelist /gen_cov/ maxtimes, do_mv, remove_mean, limit, hwindow


! ************************************************
! *** Init                                     ***
! ************************************************

  write (*,'(10x,a)') '*******************************************'
  write (*,'(10x,a)') '*             GENERATE_COVAR              *'
  write (*,'(10x,a)') '*                                         *'
  write (*,'(10x,a)') '*    Compute covariance matrix and mean   *'
  write (*,'(10x,a)') '*     state from a sequence of states.    *'
  write (*,'(10x,a)') '*                                         *'
  write (*,'(10x,a)') '*   Write covar matrix as scaled eigen-   *'
  write (*,'(10x,a)') '*    vectors and singular values into     *'
  write (*,'(10x,a)') '*              NetCDF file                *'
  write (*,'(10x,a/)') '*******************************************'

  do_mv = 0
  remove_mean = 0
  limit = 1.0e-10
  hwindow = 5

  ! Read configuration
  open(unit=20, file='pdaf.nml')
  read(unit=20, nml=gen_cov)
  rewind(20)
  close(20)

  ! Set pointer to PDAF's ensemble array
  call PDAF_set_ens_pointer(states, status)

  maxtimes = dim_ens


! ***********************************************
! *** Compute and subtract running mean state ***
! ***********************************************

   write (*,'(1x,a,i5,a)') 'Compute running mean over', maxtimes,' snapshots'
   allocate(meanstate(dim_state_p))
   allocate(tmp_ens(dim_state_p, maxtimes))
   call memcount(4,'r',dim_state_p + dim_state_p*maxtimes)

   call compute_running_mean(maxtimes, hwindow, states, meanstate, tmp_ens)

   ! get residual
   states = states - tmp_ens


   ! *********************************************************
   ! *** Singular value decomposition of covariance matrix ***
   ! ***                                                   ***
   ! *** The covariance matrix is given by the state       ***
   ! *** sequences X of k states as                        ***
   ! ***          -1    _     _ T        T                 ***
   ! *** P = (k-1)   (X-X) (X-X)  = U L U  (EVP)           ***
   ! ***                                                   ***
   ! *** Thus we compute the singular value decomposition  ***
   ! ***     _        T            -1    2  T              ***
   ! ***   X-X = U S V ;  P = (k-1)   U S  U               ***
   ! ***                                                   ***
   ! ***                         -1/2                      ***
   ! *** and we store U and (k-1)     S in a NetCDF file.  ***
   ! *********************************************************

   write (*,'(/1x,a)') '------- Compute covariance matrix decomposition -------------'

   ! PDAF_eofcovar should compute and subtract the mean state from Trajectory
   ! before decomposition. Afterwards, it's added again.


   ! Allocate arrays for singular values and vectors
   allocate(svals(maxtimes))
   allocate(stddev(n_fields))
   call memcount(4,'r',maxtimes + n_fields)


   ! Call routine generating matrix decomposition
   call PDAF_eofcovar(dim_state_p, maxtimes, n_fields, sfields(:)%dim, sfields(:)%off, &
      remove_mean, do_mv, states, stddev, svals, tmp_ens, meanstate, 1, status)

   if (status /= 0) then
      write (*, '(/1x,a)') '--------- Error in PDAF_eofcovar -----'
      stop
   end if

   ! *** determine rank to write ***
   getlimit: do i = 1, maxtimes
      if (svals(i) >= limit) then
         rank = i
      else
         exit getlimit
      end if
   end do getlimit

   if (rank < maxtimes) then
      write (*,'(1x,a,i6,a,es10.2)') &
      'use maximum of ', rank, ' eigenvectors due to eigenvalue-limit of ',limit
   end if
   if (rank == maxtimes) then
      rank = maxtimes - 1
      write (*,'(5x,a,i4)') '++ reset rank to ',rank
   end if

   write (*,'(5x,a)') 'singular values: '
   do i = 1, rank
      write (*, '(10x, i4, es12.3)') i, svals(i)
   end do


! ************************************
! *** write covariance matrix file ***
! ************************************

   call write_covar_mv(path_covar, file_covar, dim_state_p, rank, tmp_ens, meanstate, svals)


! ********************
! *** Finishing up ***
! ********************

   deallocate(meanstate)
   deallocate(tmp_ens, svals, stddev)

 contains
  

   subroutine compute_running_mean(maxtimes, hwindow, states, meanstate, run_meanstate)

     implicit none

! *** Arguments ***
     integer, intent(in) :: maxtimes         ! max number of time steps
     integer, intent(in) :: hwindow          ! half the window of running mean
     real, intent(in   ) :: states(:,:)      ! state trajectory
     real, intent(inout) :: meanstate(:)     ! mean
     real, intent(inout) :: run_meanstate(:,:) ! running means

! *** Local variables ***
     integer :: i, j                         ! counters


! *** Compute running mean about step i ***

     do i = 1, maxtimes

        meanstate = 0.0

        if (i > hwindow .and. i <= (maxtimes - hwindow)) then
           write (*,'(a,i5,1x,a,2i4)') 'step', i, 'compute mean for range', i-hwindow, i+hwindow
           do j = i - hwindow, i + hwindow
              meanstate = meanstate + states (:, j)
           end do

        else if (i <= hwindow) then
           write (*,'(a,i5,1x,a,2i4,a,2i4)') 'step', i, 'compute mean for ranges', &
                1, i+hwindow,' and ', maxtimes - hwindow +i, maxtimes
           do j = 1, i + hwindow
              meanstate = meanstate + states (:, j)
           end do
           do j = maxtimes - hwindow +i, maxtimes
              meanstate = meanstate + states (:, j)
           end do

        else if (i > maxtimes - hwindow) then
           write (*,'(a,i5,1x,a,2i4,a,2i4)') 'step', i, 'compute mean for ranges', &
                i - hwindow, maxtimes,' and ', 1, i + hwindow - maxtimes
           do j = i - hwindow, maxtimes
              meanstate = meanstate + states (:, j)
           end do
           do j = 1, i + hwindow - maxtimes
              meanstate = meanstate + states (:, j)
           end do
        end if

        meanstate = meanstate / (2.0 * real(hwindow) + 1.0)

        run_meanstate(:, i) = meanstate
     end do

   end subroutine compute_running_mean

end subroutine eofcovar
