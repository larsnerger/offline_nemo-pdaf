!$Id$
!> Building the Statevector
!!
!! This module provides variables & routines for
!! building the state vector.
!!
module mod_statevector_pdaf

  implicit none
  save

  !---- `field_ids` and `state_field` need to be adapted for a DA case -----

  ! Declare Fortran type holding the indices of model fields in the state vector
  ! This can be extended to any number of fields - it serves to give each field a name
  type field_ids
     ! Ocean Physics
     integer :: ssh
     integer :: temp
     integer :: salt
     integer :: uvel
     integer :: vvel

     ! ERGOM
     integer, allocatable  :: bgc1(:)
     integer, allocatable  :: bgc2(:)
  end type field_ids

  ! Declare Fortran type holding the definitions for model fields
  type state_field
     integer :: ndims                      ! Number of field dimensions (2 or 3)
     integer :: dim                        ! Dimension of the field
     integer :: off                        ! Offset of field in state vector
     integer :: jptrc = 0                  ! index of the tracer in nemo tracer variable
     character(len=10) :: variable = ''    ! Name of field
     character(len=20) :: name_incr = ''   ! Name of field in increment file
     character(len=20) :: name_rest_n = '' ! Name of field in restart file (n-field)
     character(len=20) :: name_rest_b = '' ! Name of field in restart file (b-field)
     character(len=30) :: file = ''        ! File name stub to read field from
     character(len=30) :: file_post = ''   ! File name part after dates
     character(len=30) :: rst_file = ''    ! Name of restart file
     character(len=20) :: unit = ''        ! Unit of variable
     integer :: transform = 0          ! Type of variable transformation
     real :: trafo_shift = 0.0         ! Constant to shift value in transformation
     integer :: limit = 0              ! Whether to limit the value of the variable
                     ! 0: no limits, 1: lower limit, 2: upper limit, 3: both limits
     real :: max_limit = 0.0           ! Upper limit of variable
     real :: min_limit = 0.0           ! Lower limit of variable
  end type state_field


  !---- The next variables usually do not need editing -----

  ! Type variable holding field IDs in state vector
  type(field_ids) :: id

  ! Type variable holding the defintions of model fields
  type(state_field), allocatable :: sfields(:)

  ! Variables to handle multiple fields in the state vector
  integer :: n_fields      !< number of fields in state vector

  ! Variables for ERGOM fields
  integer :: n_trc = 0     !< number of tracer fields
  integer :: n_bgc1 = 0    !< number of prognostic tracer fields
  integer :: n_bgc2 = 0    !< number of diagnostic tracer fields

contains

  !> This routine calculates the dimension of the local statevector.
  !!
  subroutine setup_state(dim_state_p)

    use mpi
    use mod_parallel_pdaf, &
         only: mype=>mype_ens, npes=>npes_ens, task_id, comm_ensemble, &
         comm_model, MPIerr
    use mod_nemo_pdaf, &
         only: sdim2d, sdim3d, jptra, jptra2
    use mod_memcount_pdaf, &
         only: memcount

    implicit none

! *** Arguments ***
    integer, intent(out) :: dim_state_p  !< Local dimension of state vector

! *** Local variables *** 
    integer :: i, cnt            ! Counters
    integer :: id_bgc1, id_bgc2  ! Counters
    integer :: screen=1          ! Verbosity flag
    integer :: id_var            ! Index of a variable in state vector
    integer :: dim_state         ! Global state dimension

    ! Variables to activate a field from the namelist
    ! ---- This needs to be adapted according to possible fields -----
    logical :: sv_temp = .false. ! Whether to include temperature in state vector
    logical :: sv_salt = .false. ! Whether to include salinity in state vector
    logical :: sv_ssh = .false.  ! Whether to include SSH in state vector
    logical :: sv_uvel = .false. ! Whether to include u-velocity in state vector
    logical :: sv_vvel = .false. ! Whether to include v-velocity in state vector
    logical, allocatable :: sv_bgc1(:) ! Whether to include prognostic ERGOM variables 
    logical, allocatable :: sv_bgc2(:) ! Whether to include diagnostic ERGOM variables

    ! Namelist to define active parts of state vector
    namelist /state_vector/ screen, sv_temp, sv_salt, sv_ssh, sv_uvel, sv_vvel, &
         sv_bgc1, sv_bgc2

    !---- End of section to be adapted ----


! **********************
! *** Initialization ***
! **********************

! *** Allocate arrays for ERGOM fields in state vector

    allocate(id%bgc1(jptra))
    allocate(id%bgc2(jptra2))
    id%bgc1(:)=0
    id%bgc2(:)=0

    allocate(sv_bgc1(jptra))
    allocate(sv_bgc2(jptra2))
    sv_bgc1(:) = .false.
    sv_bgc2(:) = .false.

! *** Read namelist file for state vector setup

    open (500,file='pdaf.nml')
    read (500,NML=state_vector)
    close (500)

! *** Now setup field indices in state vector

    !---- This part needs to be adapted according to possible fields in the state vector ----

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

!    if (sv_oxy) then
!       cnt = cnt + 1
!       id%oxy = cnt
!    end if

    do id_bgc1 = 1, jptra
      if (sv_bgc1(id_bgc1)) then
        cnt = cnt + 1
        id%bgc1(id_bgc1) = cnt
        n_bgc1=n_bgc1+1
      end if
    end do

    do id_bgc2 = 1, jptra2
      if (sv_bgc2(id_bgc2)) then
        cnt = cnt + 1
        id%bgc2(id_bgc2) = cnt
        n_bgc2=n_bgc2+1
      end if
    end do

    !---- End of section to be adapted ----


! ************************************************
! *** Specify state vector and state dimension ***
! ************************************************

    ! Number of model fields in state vector
    n_fields = cnt
    n_trc=n_bgc1+n_bgc2

    allocate(sfields(n_fields))


! *** Specifications for each model field in state vector ***

    !---- This part needs to be adapted according to possible fields in the state vector ----

    ! SSH
    id_var = id%ssh
    if (id_var>0) then
       sfields(id_var)%ndims = 2
       sfields(id_var)%dim = sdim2d
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

    ! Temperature
    id_var = id%temp
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
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

    ! Salinity
    id_var = id%salt
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
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

    ! U-velocity
    id_var = id%uvel
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
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

    ! V-velocity
    id_var = id%vvel
    if (id_var>0) then
       sfields(id_var)%ndims = 3
       sfields(id_var)%dim = sdim3d
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

    ! Oxygen
!     id_var = id%oxy
!     if (id_var>0) then
!        sfields(id_var)%ndims = 3
!        sfields(id_var)%dim = sdim3d
!        sfields(id_var)%variable = 'OXY'
!        sfields(id_var)%name_incr = 'bckinoxy'
!        sfields(id_var)%name_rest_n = 'TRNOXY'
!        sfields(id_var)%name_rest_b = 'TRBOXY'
!        sfields(id_var)%file = 'NORDIC_1d_ERGOM_T_'
!        sfields(id_var)%rst_file = 'restart_trc_in.nc'
!        sfields(id_var)%unit = 'mmol m-3'
!        sfields(id_var)%transform = 0
!        sfields(id_var)%trafo_shift = 0.0
!        sfields(id_var)%limit = 3
!        sfields(id_var)%min_limit = -450.0D0
!        sfields(id_var)%max_limit = 450.0D0
!     endif


    ! BGC
    do id_bgc1 = 1, jptra
      if (sv_bgc1(id_bgc1)) then
        id_var=id%bgc1(id_bgc1)
        sfields(id_var)%ndims = 3
        sfields(id_var)%dim = sdim3d
        sfields(id_var)%jptrc = id_bgc1
        sfields(id_var)%file = 'NORDIC_1d_ERGOM_T_'
        sfields(id_var)%rst_file = 'restart_trc_in.nc'
        sfields(id_var)%transform = 0
        sfields(id_var)%trafo_shift = 0.0

        select case (id_bgc1)
        case (1)
          sfields(id_var)%variable = 'NH4'
          sfields(id_var)%name_incr = 'bckinnh4'
          sfields(id_var)%name_rest_n = 'TRNNH4'
          sfields(id_var)%name_rest_b = 'TRBNH4'
          sfields(id_var)%unit = 'mmol m-3'
        case (2)
          sfields(id_var)%variable = 'NO3'
          sfields(id_var)%name_incr = 'bckinno3'
          sfields(id_var)%name_rest_n = 'TRNNO3'
          sfields(id_var)%name_rest_b = 'TRBNO3'
          sfields(id_var)%unit = 'mmol m-3'
        case (3)
          sfields(id_var)%variable = 'PO4'
          !sfields(id_var)%name_incr = 'bckinpo3'
          sfields(id_var)%name_rest_n = 'TRNPO4'
          sfields(id_var)%name_rest_b = 'TRBPO4'
          sfields(id_var)%unit = 'mmol m-3'
        case (4)
          sfields(id_var)%variable = 'SIL'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNSIL'
          sfields(id_var)%name_rest_b = 'TRBSIL'
          sfields(id_var)%unit = 'mmol m-3'
        case (5)
          sfields(id_var)%variable = 'DIA'
          sfields(id_var)%name_incr = 'bckindia'
          sfields(id_var)%name_rest_n = 'TRNDIA'
          sfields(id_var)%name_rest_b = 'TRBDIA'
          sfields(id_var)%unit = 'mmol m-3'
        case (6)
          sfields(id_var)%variable = 'FLA'
          sfields(id_var)%name_incr = 'bckinfla'
          sfields(id_var)%name_rest_n = 'TRNFLA'
          sfields(id_var)%name_rest_b = 'TRBFLA'
          sfields(id_var)%unit = 'mmol m-3'
        case (7)
          sfields(id_var)%variable = 'CYA'
          sfields(id_var)%name_incr = 'bckincya'
          sfields(id_var)%name_rest_n = 'TRNCYA'
          sfields(id_var)%name_rest_b = 'TRBCYA'
          sfields(id_var)%unit = 'mmol m-3'
        case (8)
          sfields(id_var)%variable = 'MEZ'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNMEZ'
          sfields(id_var)%name_rest_b = 'TRBMEZ'
          sfields(id_var)%unit = 'mmol m-3'
        case (9)
          sfields(id_var)%variable = 'MIZ'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNMIZ'
          sfields(id_var)%name_rest_b = 'TRBMIZ'
          sfields(id_var)%unit = 'mmol m-3'
        case (10)
          sfields(id_var)%variable = 'DET'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNDET'
          sfields(id_var)%name_rest_b = 'TRBDET'
          sfields(id_var)%unit = 'mmol m-3'
        case (11)
          sfields(id_var)%variable = 'DETs'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNDETs'
          sfields(id_var)%name_rest_b = 'TRBDETs'
          sfields(id_var)%unit = 'mmol m-3'
        case (12)
          sfields(id_var)%variable = 'FE'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNFE'
          sfields(id_var)%name_rest_b = 'TRBFE'
          sfields(id_var)%unit = 'mmol m-3'
        case (13)
          sfields(id_var)%variable = 'LDON'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNLDON'
          sfields(id_var)%name_rest_b = 'TRBLDON'
          sfields(id_var)%unit = 'mmol m-3'
        case (14)
          sfields(id_var)%variable = 'DIC'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNDIC'
          sfields(id_var)%name_rest_b = 'TRBDIC'
          sfields(id_var)%unit = 'mmol m-3'
        case (15)
          sfields(id_var)%variable = 'ALK'
          !sfields(id_var)%name_incr = ''
          sfields(id_var)%name_rest_n = 'TRNALK'
          sfields(id_var)%name_rest_b = 'TRBALK'
          sfields(id_var)%unit = 'mmol m-3'
        case (16)
          sfields(id_var)%variable = 'OXY'
          sfields(id_var)%name_incr = 'bckinoxy'
          sfields(id_var)%name_rest_n = 'TRNOXY'
          sfields(id_var)%name_rest_b = 'TRBOXY'
          sfields(id_var)%unit = 'mmol m-3'
        end select
      end if
    end do

    do id_bgc2 = 1, jptra2
      if (sv_bgc2(id_bgc2)) then
        id_var=id%bgc2(id_bgc2)
        sfields(id_var)%ndims = 3
        sfields(id_var)%dim = sdim3d
        sfields(id_var)%jptrc = id_bgc2
        sfields(id_var)%file = 'NORDIC_1d_ERGOM_T_'
        !sfields(id_var)%rst_file = ''
        sfields(id_var)%transform = 0
        sfields(id_var)%trafo_shift = 0.0

        select case (id_bgc2)
        case (1)
          sfields(id_var)%variable = 'xpco2'
          !sfields(id_var)%name_incr = ''
          !sfields(id_var)%name_rest_n = ''
          !sfields(id_var)%name_rest_b = ''
          sfields(id_var)%unit = 'micro atm'
        case (2)
          sfields(id_var)%variable = 'xph'
          !sfields(id_var)%name_incr = ''
          !sfields(id_var)%name_rest_n = ''
          !sfields(id_var)%name_rest_b = ''
          sfields(id_var)%unit = '-'
        case (3)
          sfields(id_var)%variable = 'xchl'
          !sfields(id_var)%name_incr = ''
          !sfields(id_var)%name_rest_n = ''
          !sfields(id_var)%name_rest_b = ''
          sfields(id_var)%unit = 'mg m-3'
        end select
      end if
    end do

    !---- End of section to be adapted ----


! *** Compute offsets ***

    ! Define offsets in state vector
    sfields(1)%off = 0
    do i = 2, n_fields
       sfields(i)%off = sfields(i-1)%off + sfields(i-1)%dim
    end do

! *** Set state vector dimension ***

    dim_state_p = sum(sfields(:)%dim)


! *** Write information about the state vector ***

    if (mype==0) then 
       write (*,'(/a,2x,a)') 'NEMO-PDAF', '*** Setup of state vector ***' 
       write (*,'(a,5x,a,i5)') 'NEMO-PDAF', '--- Number of fields in state vector:', n_fields
       if (npes==1) then
          write (*,'(a,5x,a2,3x,a8,6x,a5,7x,a3,7x,a6)') 'NEMO-PDAF','ID', 'variable', 'ndims', 'dim', 'offset'
          write (*,'(a,5x,49a)') 'NEMO-PDAF', ('-',i=1,49)
          do i = 1, n_fields
             write (*,'(a, 2x,i5,3x,a10,2x,i5,3x,i10,3x,i10)') &
                  'NEMO-PDAF',i, sfields(i)%variable, sfields(i)%ndims, sfields(i)%dim, sfields(i)%off
          end do
       else
          if (task_id==1) then
             write (*,'(a,a4,7x,a2,3x,a8,6x,a5,7x,a3,7x,a6)') &
                  'NEMO-PDAF','pe','ID', 'variable', 'ndims', 'dim', 'offset'
          end if
       end if
    end if

    if (npes>1 .and. (mype==0 .or. (task_id==1 .and. screen>2))) then
       do i = 1, n_fields
          write (*,'(a,i3,5x,i5,3x,a10,2x,i5,3x,i10,3x,i10)') &
               'NEMO-PDAF',mype, i, sfields(i)%variable, sfields(i)%ndims, sfields(i)%dim, sfields(i)%off
       end do
    end if

    if (npes==1) then
       write (*,'(a,2x,a,1x,i10)') 'NEMO-PDAF', 'Full state dimension: ',dim_state_p
    else
       if (task_id==1) then
          if (screen>1 .or. mype==0) &
               write (*,'(a,2x,a,1x,i4,2x,a,1x,i10)') &
               'NEMO-PDAF', 'PE', mype, 'PE-local full state dimension: ',dim_state_p

          call MPI_Reduce(dim_state_p, dim_state, 1, MPI_INTEGER, MPI_SUM, 0, COMM_model, MPIerr)
          if (mype==0) then
             write (*,'(a,2x,a,1x,i10)') 'NEMO-PDAF', 'Global state dimension: ',dim_state
          end if
       end if
    end if
    call MPI_Barrier(comm_ensemble, MPIerr)


  end subroutine setup_state

end module mod_statevector_pdaf
