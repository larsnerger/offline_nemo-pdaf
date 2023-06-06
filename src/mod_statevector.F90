!$Id$
!BOP
!
! !MODULE:
module mod_statevector

! !DESCRIPTION:
! This module provides variables for the setup of
! the state vector.
!
! !REVISION HISTORY:
! 2022-04 - Lars Nerger - Initial code by splitting mod_assimilation
! Later revisions - see svn log
!
! !USES:
  implicit none
  save
!EOP

  integer :: n_fields      !< number of fields in state vector

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
     integer :: oxy
  end type field_ids

  ! Type variable holding field IDs in state vector
  type(field_ids) :: id

  type state_field
     integer :: ndims                  ! Number of field dimensions (2 or 3)
     integer :: dim                    ! Dimension of the field
     integer :: off                    ! Offset of field in state vector
     character(len=10) :: variable     ! Name of field
     character(len=20) :: name_incr    ! Name of field in increment file
     character(len=20) :: name_rest_n  ! Name of field in restart file (n-field)
     character(len=20) :: name_rest_b  ! Name of field in restart file (b-field)
     character(len=30) :: file         ! File name stub to read field from
     character(len=30) :: file_post='' ! File name part after dates
     character(len=30) :: rst_file     ! Name of restart file
     character(len=20) :: unit         ! Unit of variable
     integer :: transform = 0          ! Type of variable transformation
     real :: trafoConst = 0.0          ! Constant to shift value in transformation
     integer :: limit = 0              ! Whether to limit the value of the variable
                                       ! 0: no limits, 1: lower limit, 2: upper limit, 3: both limits
     real :: max_limit = 0.0           ! Upper limit of variable
     real :: min_limit = 0.0           ! Lower limit of variable
  end type state_field

  type(state_field), allocatable :: sfields(:)

end module mod_statevector
