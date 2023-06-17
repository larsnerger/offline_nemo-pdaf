module mod_aux_pdaf

! !DESCRIPTION:
! Auxiliary routines
! - mapping between state vector and model field
! - state transformations (log, etc.)

  ! Include dimension information for model grid
  use mod_nemo_pdaf, &
       only: nlvls=>jpk, nj_p, ni_p, nwet, wet_pts, use_wet_state
  use mod_statevector_pdaf, &
       only: id, sfields, n_fields

  implicit none
  save
  
contains
!===============================================================================


!> Convert from NEMO model field to state vector
!!
  subroutine field2state(field, state, offset, ndims, missval)

    implicit none

! *** Arguments ***
    real, intent(in)    :: field(:,:,:,:)   !< Model field
    real, intent(inout) :: state(:)         !< State vector
    integer, intent(in) :: offset           !< Offset in state vector
    integer, intent(in) :: ndims            !< Number of dimensions in field
    real, intent(in) :: missval             !< missing value

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize state vector from model field ***

    if (use_wet_state==1) then
       
       do k = 1, n_levels
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             if (k <= wet_pts(3,i)) then
                state(cnt) = field(wet_pts(6, i), wet_pts(7, i), k, 1)
             end if
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then

!$OMP PARALLEL DO PRIVATE (i, k, cnt)
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                state(cnt + k) = field(wet_pts(6, i), wet_pts(7, i), k, 1)
             end do
          end do
!$OMP END PARALLEL DO

       else

!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = offset + i
             state(cnt) = field(wet_pts(6, i), wet_pts(7, i), 1, 1)
          end do
!$OMP END PARALLEL DO

       end if

    else
       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1,nj_p
             do i = 1, ni_p
                if (abs(field(i,j,k,1)- missval) > 0.1) then 
                   state(cnt) = field(i,j,k,1)
                else
                   state(cnt) = 0.0
                endif
                cnt = cnt + 1
             enddo
          enddo
       enddo
    end if

  end subroutine field2state

!==============================================================================

!> Convert from state vector to NEMO model field
!!
  subroutine state2field(state, field, offset, ndims)

    implicit none

! *** Arguments ***
    real, intent(in)    :: state(:)         !< State vector
    real, intent(out)   :: field(:,:,:,:)   !< Model field
    integer, intent(in) :: offset           !< Offset in state vector
    integer, intent(in) :: ndims            !< Number of dimensions in field

! *** Local variables ***
    integer :: i, j, k
    integer :: cnt
    integer :: n_levels


! *** Set number of model layers ***
    if (ndims == 3) then
       n_levels = nlvls
    else
       n_levels = 1
    end if


! *** Initialize model field from state vector

    if (use_wet_state==1) then
       
       do k = 1, nlvls
!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = i + nwet*(k-1) + offset
             field(wet_pts(6, i), wet_pts(7, i), k, 1) = state(cnt)
          end do
!$OMP END PARALLEL DO
       end do

    elseif (use_wet_state==2) then

       if (ndims == 3) then

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = wet_pts(5,i) - 1 + offset
             do k = 1, wet_pts(3,i)
                field(wet_pts(6, i), wet_pts(7, i), k, 1) = state(cnt + k)
             end do
!!!!$OMP END PARALLEL DO
          end do

       else

!!!!$OMP PARALLEL DO PRIVATE (i, cnt)
          do i = 1, nwet
             cnt = offset + i
             field(wet_pts(6, i), wet_pts(7, i), 1, 1) = state(cnt)
!!!!$OMP END PARALLEL DO
          end do

       end if
    else

       cnt = 1 + offset
       do k = 1, n_levels
          do j = 1, nj_p
             do i = 1, ni_p
                field(i,j,k,1) = state(cnt) !convert to NEMO format (ntimec,nlvls,nlats,nlons)
                cnt = cnt + 1
             enddo
          enddo
       enddo

    end if

  end subroutine state2field

!==============================================================================

!> Transform field, e.g. to log and back
!!
  subroutine transform_field(type, trafo, shift, state, dim, off, var)

    implicit none

    integer, intent(in) :: type     !< Direction transformation
    integer, intent(in) :: trafo    !< Type of transformation
    real, intent(in)    :: shift    !< constant for shifting value in transformation
    real, intent(inout) :: state(:) !< State vector
    integer, intent(in) :: dim      !< dimension of field in state vector
    integer, intent(in) :: off      !< Offset of field in state vector
    character(len=*),intent(in) :: var      !< Name of variable


    if (type==1) then
       ! Transformation from NEMO value to transformed value

       select case (trafo)
       case(0)
!          write(*,*) 'No Transformation of variable ', trim(var)
       case(1)
          write(*,*) 'use log basis 10 transformation of variable ', trim(var)
          state(off+1 : off+dim) = log10(state(off+1 : off+dim) + shift)
       case(2)
          write(*,*) 'use ln transformation of variable ', trim(var)
          state(off+1 : off+dim) = log(state(off+1 : off+dim) + shift)
       case(3)
          write(*,*) 'no transformation- box cox still needs to be implemented'
       case DEFAULT
!          write(*,*) 'No Transformation of variable ', trim(var)
       end select

    elseif (type==2) then
       ! Transformation from transformed value to NEMO value

       select case (trafo)
       case(0)
!          write(*,*) 'No Transformation of bio limit'
       case(1)
          write(*,*) 'use log basis 10 transformation of bio variable'
          state(off+1 : off+dim) = (10.D0**state(off+1 : off+dim))-shift
       case(2)
          write(*,*) 'use ln transformation of bio variable'
          state(off+1 : off+dim) = (exp(state(off+1 : off+dim)))-shift
       case(3)
          write(*,*) 'no transformation- box cox still needs to be implemented'
       case DEFAULT
!          write(*,*) 'No Transformation of bio variable'
       end select

    end if

  end subroutine transform_field

!==============================================================================

!> Transform all fields, e.g. to log and back
!!
  subroutine transform_field_mv(type, state)

    implicit none

! *** Arguments ***
    integer, intent(in) :: type     !< Direction of transformation
    real, intent(inout) :: state(:) !< State vector

! *** Local variables ***
    integer :: i, j       ! Counters
    integer :: trafo      ! Type of transformation
    real    :: shift      ! constant for shifting value in transformation
    integer :: dim        ! dimension of field in state vector
    integer :: off        ! Offset of field in state vector
    character(len=10) :: var      ! Name of variable
    integer :: dolimit    ! Whether to apply a min/max limit
    real    :: max_limit  ! Maximum limiting value
    real    :: min_limit  ! Minimum limiting value


    do i = 1, n_fields

       ! Initialize values from sfields
       trafo = sfields(i)%transform
       shift = sfields(i)%trafo_shift
       dim = sfields(i)%dim
       off = sfields(i)%off
       var = sfields(i)%variable
       dolimit = sfields(i)%limit
       min_limit = sfields(i)%min_limit
       max_limit = sfields(i)%max_limit

       
       ! *** Apply Limits
       if (type==2) then
          call var_limits_mv(state)
       end if

       ! *** Transformations

       if (type==1) then
          ! Transformation from NEMO value to transformed value

          select case (trafo)
          case(0)
!          write(*,*) 'No Transformation of variable ', trim(var)
          case(1)
             write(*,*) 'use log basis 10 transformation of variable ', trim(var)
             state(off+1 : off+dim) = log10(state(off+1 : off+dim) + shift)
          case(2)
             write(*,*) 'use ln transformation of variable ', trim(var)
             state(off+1 : off+dim) = log(state(off+1 : off+dim) + shift)
          case(3)
             write(*,*) 'no transformation- box cox still needs to be implemented'
          case DEFAULT
!          write(*,*) 'No Transformation of variable ', trim(var)
          end select

       elseif (type==2) then
          ! Transformation from transformed value to NEMO value

          select case (trafo)
          case(0)
!          write(*,*) 'No Transformation of bio limit'
          case(1)
             write(*,*) 'use log basis 10 transformation of variable ', trim(var)
             state(off+1 : off+dim) = (10.D0**state(off+1 : off+dim))-shift
          case(2)
             write(*,*) 'use ln transformation of variable ', trim(var)
             state(off+1 : off+dim) = (exp(state(off+1 : off+dim)))-shift
          case(3)
             write(*,*) 'no transformation- box cox still needs to be implemented'
          case DEFAULT
!          write(*,*) 'No Transformation of bio variable'
          end select

       end if

    end do

  end subroutine transform_field_mv

!==============================================================================

!> Apply limits to variables
!!
  subroutine var_limits_mv(state)

    implicit none

! *** Arguments ***
    real, intent(inout) :: state(:) !< State vector

! *** Local variables ***
    integer :: i, j       ! Counters
    integer :: dim        ! dimension of field in state vector
    integer :: off        ! Offset of field in state vector
    character(len=10) :: var      ! Name of variable
    integer :: dolimit    ! Whether to apply a min/max limit
    real    :: max_limit  ! Maximum limiting value
    real    :: min_limit  ! Minimum limiting value


    do i = 1, n_fields

       ! Initialize values from sfields
       dim = sfields(i)%dim
       off = sfields(i)%off
       var = sfields(i)%variable
       dolimit = sfields(i)%limit
       min_limit = sfields(i)%min_limit
       max_limit = sfields(i)%max_limit


       ! *** Apply Limits

       if (dolimit == 1) then

          ! Apply minimum limit
          do j = off+1, off+dim
             if (state(j) < min_limit) state(j) = min_limit
          end do

       elseif (dolimit == 2) then

          ! Apply maximum limit
          do j = off+1, off+dim
             if (state(j) > max_limit) state(j) = max_limit
          end do

       elseif (dolimit == 3) then

          ! Apply minimum and maximum limits
          do j = off+1, off+dim
             if (state(j) < min_limit) then
                state(j) = min_limit
             elseif (state(j) > max_limit) then
                state(j) = max_limit
             end if
          end do

       end if

    end do

  end subroutine var_limits_mv

end module mod_aux_pdaf
