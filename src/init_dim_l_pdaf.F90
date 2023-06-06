!$Id$
!>  Set dimension of local model state
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during analysis step
!! in PDAF_X_update in the loop over all local
!! analysis domains. It has to set the dimension
!! of the local model  state on the current analysis
!! domain.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
subroutine init_dim_l_pdaf(step, domain_p, dim_l)

   use mod_assimilation_pdaf, &
        only:  domain_coords, id_lstate_in_pstate, dim_state_p
   use mod_statevector_pdaf, &
        only: n_fields, id, sfields
   use mod_nemo_pdaf, &
        only:  nwet, wet_pts, dim_2d, nlons=>jpiglo, nlats=>jpjglo, &
        use_wet_state

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step     !< Current time step
  integer, intent(in)  :: domain_p !< Current local analysis domain
  integer, intent(out) :: dim_l    !< Local state dimension

! *** Local variables ***
  integer(4) :: i, cnt, ifield
  integer(4) :: domain_all
  integer(4) :: id_surf

  
! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  ! Distinguish 2D and 3D fields
  dim_l = 0
  do i = 1, n_fields
     if (sfields(i)%ndims==2) then
        dim_l = dim_l + 1
     else
        dim_l = dim_l + wet_pts(3, domain_p)
     end if
  end do


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Initialize domain index in nlons*nlats domain
  domain_all = wet_pts(4, domain_p)

  ! Set coordinates according to grid point indices
  domain_coords(1) = real(wet_pts(1, domain_p))
  domain_coords(2) = real(wet_pts(2, domain_p))


! ******************************************************
! *** Initialize array of indices of the local state ***
! ***  vector elements in the global state vector.   ***
! ******************************************************

  ! Allocate array
  if (allocated(id_lstate_in_pstate)) deallocate(id_lstate_in_pstate)
  allocate(id_lstate_in_pstate(dim_l))

  cnt = 1

  do ifield = 1, n_fields

     ! Set indices
     if (use_wet_state==1) then
        ! state vector contains full columns of surface wet points

        id_surf = domain_p + sfields(ifield)%off

        id_lstate_in_pstate(cnt) = id_surf
        cnt = cnt + 1

        do i = 2, wet_pts(3,domain_p)
           id_lstate_in_pstate(cnt) = id_surf + (i-1)*nwet
           cnt = cnt + 1
        enddo

     elseif (use_wet_state==2) then
        ! state vector only contains wet points - stored in leading vertical order

        if (sfields(ifield)%ndims==3) then

           ! 3D field
           id_surf = wet_pts(5, domain_p) + sfields(ifield)%off

           do i = 1, wet_pts(3,domain_p)
              id_lstate_in_pstate(cnt) = id_surf + i - 1
              cnt = cnt  + 1
           enddo

        else

           ! 2D field
           id_lstate_in_pstate(cnt) = domain_p + sfields(ifield)%off
           cnt = cnt  + 1

        end if

     else
        ! State vector contains full 3d box

        ! This needs to be checked for MPI -decomposition!

        id_surf = mod(domain_all, dim_2d) + sfields(ifield)%off

        id_lstate_in_pstate(cnt) = id_surf
        cnt = cnt + 1
        do i = 2, wet_pts(3,domain_p)
           id_lstate_in_pstate(cnt) = id_surf + (i-1)*dim_2d
           cnt = cnt + 1
        enddo

     end if
  
  enddo

  if (id_lstate_in_pstate(wet_pts(3,domain_p)) > dim_state_p) then
     write(*,*) 'Error: please check the global indices for local state vector'
  endif

end subroutine init_dim_l_pdaf
