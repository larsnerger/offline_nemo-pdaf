!$Id$
!> Initialize setup of state vector
!!
!! Routine to perform initialize the multivarate state vector.
!! Here the different fields contained in the state vector are defined.
!!
!! !REVISION HISTORY:
!! 2022-02 - Lars Nerger - Initial code
!! Later revisions - see svn log
!!
subroutine setup_state(dim_state_p)

  use mpi
  use mod_parallel_pdaf, &
       only: mype=>mype_world, npes=>npes_world, MPIerr
  use mod_statevector_pdaf, &
       only: n_fields, sfields, id
  use mod_nemo_pdaf, &
       only: jpiglo, jpjglo, jpk, nwet, nwet3d, use_wet_state
  use mod_memcount_pdaf, &
       only: memcount

  implicit none

! *** Arguments ***
  integer, intent(out) :: dim_state_p

! *** local variables *** 
  integer :: i, cnt            ! Counters
  integer :: id_var            ! Index of a variable in state vector
  integer :: dim2d, dim3d      ! Dimension of 2D and 3D field in state vector
  logical :: sv_temp = .false. ! Whether to include temperature in state vector
  logical :: sv_salt = .false. ! Whether to include salinity in state vector
  logical :: sv_ssh = .false.  ! Whether to include SSH in state vector
  logical :: sv_uvel = .false. ! Whether to include u-velocity in state vector
  logical :: sv_vvel = .false. ! Whether to include v-velocity in state vector
  logical :: sv_oxy = .false.  ! Whether to include oxygen in state vector

  namelist /state_vector/ sv_temp, sv_salt, sv_ssh, sv_uvel, sv_vvel, &
       sv_oxy


! **********************
! *** Initialization ***
! **********************

! *** Read namelist file for state vector setup

  open (500,file='pdaf.nml')
  read (500,NML=state_vector)
  close (500)

! *** Now setup state vector

  cnt = 0
  if (sv_ssh) then
     cnt = cnt + 1
     id%ssh = cnt
  end if

  if (sv_temp) then
     cnt = cnt + 1
     id%temp = cnt
  end if

  if (sv_salt) then
     cnt = cnt + 1
     id%salt = cnt
  end if

  if (sv_uvel) then
     cnt = cnt + 1
     id%uvel = cnt
  end if

  if (sv_vvel) then
     cnt = cnt + 1
     id%vvel = cnt
  end if

!  if (sv_oxy) then
!     cnt = cnt + 1
!     id%oxy = cnt
!  end if



! ************************************************
! *** Specify state vector and state dimension ***
! ************************************************

  ! Number of model fields in state vector
  n_fields = cnt

! *** set dimension of 2d and 3d fields in state vector ***
  
  if (use_wet_state==1) then
     ! State vector contains full columns when surface grid point is wet
     dim3d = nwet*jpk
     dim2d = nwet
  elseif (use_wet_state==2) then
     ! State vector only contains wet grid points
     dim3d = nwet3d
     dim2d = nwet
  else
     ! State vector contains 3d grid box
     dim3d = jpiglo*jpjglo*jpk
     dim2d = jpiglo*jpjglo
  end if

  allocate(sfields(n_fields))


! *** Specifications for each model field in state vector ***

  id_var = id%ssh
  if (id_var>0) then
     sfields(id_var)%ndims = 2
     sfields(id_var)%dim = dim2d
     sfields(id_var)%variable = 'SSH_inst'
     sfields(id_var)%name_incr = 'bckineta'
     sfields(id_var)%name_rest_n = 'sshn'
     sfields(id_var)%name_rest_b = 'sshb'
     sfields(id_var)%file = 'NORDIC_1d_SURF_grid_T_'
     sfields(id_var)%rst_file = 'restart_in.nc'
     sfields(id_var)%unit = 'm'
     sfields(id_var)%transform = 0
     sfields(id_var)%trafo_shift = 0.0
     sfields(id_var)%limit = 0
  endif

  id_var = id%temp
  if (id_var>0) then
     sfields(id_var)%ndims = 3
     sfields(id_var)%dim = dim3d
     sfields(id_var)%variable = 'votemper'
     sfields(id_var)%name_incr = 'bckint'
     sfields(id_var)%name_rest_n = 'tn'
     sfields(id_var)%name_rest_b = 'tb'
     sfields(id_var)%file = 'NORDIC_1d_grid_T_'
     sfields(id_var)%rst_file = 'restart_in.nc'
     sfields(id_var)%unit = 'degC'
     sfields(id_var)%transform = 0
     sfields(id_var)%trafo_shift = 0.0
  endif

  id_var = id%salt
  if (id_var>0) then
     sfields(id_var)%ndims = 3
     sfields(id_var)%dim = dim3d
     sfields(id_var)%variable = 'vosaline'
     sfields(id_var)%name_incr = 'bckins'
     sfields(id_var)%name_rest_n = 'sn'
     sfields(id_var)%name_rest_b = 'sb'
     sfields(id_var)%file = 'NORDIC_1d_grid_T_'
     sfields(id_var)%rst_file = 'restart_in.nc'
     sfields(id_var)%unit = '1e-3'
     sfields(id_var)%transform = 0
     sfields(id_var)%trafo_shift = 0.0
  endif

  id_var = id%uvel
  if (id_var>0) then
     sfields(id_var)%ndims = 3
     sfields(id_var)%dim = dim3d
     sfields(id_var)%variable = 'uos'
     sfields(id_var)%name_incr = 'bckinu'
     sfields(id_var)%name_rest_n = 'un'
     sfields(id_var)%name_rest_b = 'ub'
     sfields(id_var)%file = 'NORDIC_1d_grid_U_'
     sfields(id_var)%rst_file = 'restart_in.nc'
     sfields(id_var)%unit = 'm/s'
     sfields(id_var)%transform = 0
     sfields(id_var)%trafo_shift = 0.0
  endif

  id_var = id%vvel
  if (id_var>0) then
     sfields(id_var)%ndims = 3
     sfields(id_var)%dim = dim3d
     sfields(id_var)%variable = 'vos'
     sfields(id_var)%name_incr = 'bckinv'
     sfields(id_var)%name_rest_n = 'vn'
     sfields(id_var)%name_rest_b = 'vb'
     sfields(id_var)%file = 'NORDIC_1d_grid_V_'
     sfields(id_var)%rst_file = 'restart_in.nc'
     sfields(id_var)%unit = 'm/s'
     sfields(id_var)%transform = 0
     sfields(id_var)%trafo_shift = 0.0
  endif

!  id_var = id%oxy
  id_var = 0
  if (id_var>0) then
     sfields(id_var)%ndims = 3
     sfields(id_var)%dim = dim3d
     sfields(id_var)%variable = 'OXY'
     sfields(id_var)%name_incr = 'bckinoxy'
     sfields(id_var)%name_rest_n = 'TRNOXY'
     sfields(id_var)%name_rest_b = 'TRBOXY'
     sfields(id_var)%file = 'NORDIC_1d_ERGOM_T_'
     sfields(id_var)%rst_file = 'restart_trc_in.nc'
     sfields(id_var)%unit = 'mmol m-3'
     sfields(id_var)%transform = 0
     sfields(id_var)%trafo_shift = 0.0
     sfields(id_var)%limit = 3
     sfields(id_var)%min_limit = -450.0D0
     sfields(id_var)%max_limit = 450.0D0
  endif

! *** Compute offsets ***

  ! Define offsets in state vector
  sfields(1)%off = 0
  DO i = 2, n_fields
     sfields(i)%off = sfields(i-1)%off + sfields(i-1)%dim
  END DO

! *** Set state vector dimension ***

  dim_state_p = sum(sfields(:)%dim)


! *** Write information about the state vector ***

  if (mype==0) then 
     write (*,'(/5x,a)') '*** Setup of state vector ***' 
     write (*,'(5x,a,i5)') '--- Number of fields in state vector:', n_fields
     if (npes==1) THEN
        write (*,'(8x,a2,3x,a8,6x,a5,7x,a3,7x,a6)') 'ID', 'variable', 'ndims', 'dim', 'offset'
        write (*,'(8x,49a)') ('-',i=1,49)
        do i = 1, n_fields
           write (*,'(5x,i5,3x,a10,2x,i5,3x,i10,3x,i10)') &
                i, sfields(i)%variable, sfields(i)%ndims, sfields(i)%dim, sfields(i)%off
        end do
     else
        write (*,'(a4,7x,a2,3x,a8,6x,a5,7x,a3,7x,a6)') 'pe','ID', 'variable', 'ndims', 'dim', 'offset'
     end if
  end if
  if (npes>1) then
!     write (*,'(i3,5x,a,i12)') mype,'--- state vector dimension', dim_state_p
     do i = 1, n_fields
        write (*,'(i3,5x,i5,3x,a10,2x,i5,3x,i10,3x,i10)') &
             mype, i, sfields(i)%variable, sfields(i)%ndims, sfields(i)%dim, sfields(i)%off
     end do
  end if

  if (npes==1) then
     write (*,'(5x,a,1x,i10)') 'Full state dimension: ',dim_state_p
  else
     write (*,'(2x,a,1x,i4,2x,a,1x,i10)') 'PE', mype, 'Full state dimension: ',dim_state_p
  end if

  CALL MPI_Barrier(MPI_COMM_WORLD, MPIerr)

end subroutine setup_state
