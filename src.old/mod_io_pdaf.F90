!> Module holding IO operations for NEMO-PDAF
!!
module mod_io_pdaf

  use mpi

  ! Include dimension information for model grid
  use mod_nemo_pdaf, &
       only: nlvls=>jpk, nlats=>jpjglo, nlons=>jpiglo, &
       depths=>gdept_1d, lons=>glamt, lats=>gphit, &
       tmp_4d, ni_p, nj_p, nk_p, istart, jstart

  ! Include information on state vector
  use mod_statevector_pdaf, &
       only: id, sfields, n_fields

  ! Include parallelization information
  use mod_parallel_pdaf, &
       only: mype=>mype_ens, npes=>npes_ens, comm_filter

  ! Include auxiliary routines
  use mod_aux_pdaf, &
       only: field2state, state2field, transform_field, transform_field_mv

  implicit none
  save

  integer :: verbose=1   ! Set verbosity of IO routines (0,1,2,3)

  ! Control of IO
  character(len=4) :: write_var_time='fcst'  ! Write variance at 'fcst', 'ana', or 'none'
  logical :: write_ens_states=.false.        ! Write a singel file of ensmeble state vectors
  logical :: write_ens_fields=.false.        ! Write set of files holding ensemble fields
  logical :: saveState                       ! Write analysis state to file
  logical :: saveIncr                        ! Write increment to file
  logical :: do_deflate=.false.              ! Deflate variables in NC files (this seems to fail for parallel nc)
  character(len=4)   :: coupling_nemo = 'rest'   ! offline: 'rest', 'incr', online: 'oinc', 'odir'

  character(len=100) :: file_PDAF_state    ! File name for outputs of analysis state as fields
  character(len=100) :: file_PDAF_incr     ! File name for increment
  character(len=100) :: file_PDAF_variance ! File name for variance
  character(len=200) :: path_state         ! Path to NEMO files
  character(len=8)   :: file_state_date1   ! Date 1 in NEMO file name
  character(len=8)   :: file_state_date2   ! Date 2 in NEMO file name
  character(len=200) :: path_ens           ! Path of ensemble file  
  character(len=80)  :: file_ens           ! File name of ensemble file
  character(len=200) :: path_covar         ! Path of file holding covariance matrix
  character(len=80)  :: file_covar         ! Filename for covariance matrix
  character(len=200) :: path_restart       ! Path of restart file
  character(len=80)  :: file_restart       ! file name of restart dile
  character(len=80)  :: ens_filelist       ! Name of file holding dates to read in ensemble states
  character(len=1)   :: datestype='-'      ! Linking character in dates (`_` or `-`)

  real :: startEnsTime=1.0, endEnsTime=1.0, incrTime=1.0

  ! NEMO output file
  integer(4)        :: ntimec=1
  integer(4)        :: nfiles

  ! Missing value in netcdf file
  real(8) :: missing_value
    
contains
! ===================================================================================

!> Read fields from NEMO file into a state vector
!!
  subroutine read_state_mv(path, date1, date2, dim_p, itime, coupling, state)

    use netcdf
    
    implicit none
    
    character(len = *), intent(in)    :: path          !< Path to file
    character(len = *), intent(in)    :: date1         !< Start date in file name
    character(len = *), intent(in)    :: date2         !< End date in file name
    integer(4),         intent(in)    :: dim_p         !< PE-local state dimension
    integer(4),         intent(in)    :: itime         !< Time to read in file
    character(len = *), intent(in)    :: coupling      !< Type of NEMO coupling
    real(8),            intent(inout) :: state(dim_p)  !< State vector
    
    ! Local variables 
    integer(4) :: i, cnt          ! Counters
    integer(4) :: varid           ! Variable ID
    integer(4) :: ncid            ! NC file id
    character(len=50) :: filename ! Full file name
    character(len=17) :: dates    ! String with initial and final date

    if (verbose>0 .and. mype==0) &
         write(*,'(a,4x,a,i8)') 'NEMO-PDAF', '*** Read model output at time step: ', itime

    if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

    ! Initialize state
    state = 0.0

    if (datestype=='_') then
       dates = trim(date1)//'_'//trim(date2)
    else
       dates = trim(date1)//'-'//trim(date2)
    end if

    do i = 1, n_fields
    
       filename = trim(sfields(i)%file)//trim(dates)//trim(sfields(i)%file_post)//'.nc'
       if (verbose>1 .and. mype==0) then 
          write(*,'(a,2x,a)') 'NEMO-PDAF', trim(path)//trim(filename)
          write (*,'(a,i5,a,a,a,i10)') &
               'NEMO-PDAF', i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       ! Open the file
       call check( nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid) )

       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )
       if (coupling/='rest') then
          call check( nf90_get_att(ncid, varid, 'missing_value', missing_value) )
       else
          missing_value=0.0
       endif

       ! Read variable
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )

       if (sfields(i)%ndims == 3) then
          call check( nf90_get_var(ncid, varid, tmp_4d, &
               start=(/istart, jstart, 1, itime/), count=(/ni_p, nj_p, nlvls, 1/)) )
       else
          call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
               start=(/istart, jstart, itime/), count=(/ni_p, nj_p, 1/)) )
       end if

       call check( nf90_close(ncid) )

       ! Convert field to state vector
       call field2state(tmp_4d, state, sfields(i)%off, sfields(i)%ndims, missing_value)

    end do
   
    ! Potentially transform fields
    call transform_field_mv(1, state)

    if (verbose>2) then
       do i = 1, n_fields
          write(*,'(a, 1x, a, a10, 1x, a,5x, 2f12.6)') &
               'NEMO-PDAF', 'Min and max for ',trim(sfields(i)%variable),' :     ',              &
               minval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim)), &
               maxval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim))
       enddo
    end if

  end subroutine read_state_mv


!=============================================================================== 

!> Read an ensemble of model fields into the ensemble array
!!  
  subroutine read_ens_mv_loop(path, date1, date2, dim_p, dim_ens, coupling, ens)

    use netcdf
    
    implicit none

! *** Arguments ***    
    character(len = *), intent(in)   :: path      !< Path of file
    character(len = *), intent(in)   :: date1     !< Start date in file name
    character(len = *), intent(in)   :: date2     !< End date in file name
    integer(4),         intent(in)   :: dim_p     !< State dimension
    integer(4),         intent(in)   :: dim_ens   !< Ensemble size
    character(len = *), intent(in)   :: coupling  !< Type of NEMO coupling
    real(8),            intent(inout):: ens(:,:)  !< Ensemble array
    
! *** Local variables ***
    integer(4) :: i, cnt, member   ! Counters
    integer(4) :: ncid             ! NC file ID
    integer(4) :: varid            ! Variable ID
    character(len=50) :: filename  ! Full file name
    character(len=17) :: dates     ! Combined date string of file

    if (verbose>0 .and. mype==0) &
         write(*,'(a,4x,a)') 'NEMO-PDAF','*** Ensemble: Read model snapshots'

    if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

    ! Initialize ensemble
    ens = 0.0

    do i = 1, n_fields

       if (datestype=='_') then
          dates = trim(date1)//'_'//trim(date2)
       else
          dates = trim(date1)//'-'//trim(date2)
       end if
       filename = trim(sfields(i)%file)//trim(dates)//trim(sfields(i)%file_post)//'.nc'
       if (verbose>1 .and. mype==0) then
          write(*,'(a,2x,a)') 'NEMO-PDAF', trim(path)//trim(filename)
          write (*,'(a,i5,a,a,a,i10)') &
               'NEMO-PDAF', i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       ! Open the file
       call check( nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid) )

       !  Read field
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )

       ! Read missing value
       if (coupling/='rest') then
          call check( nf90_get_att(ncid, varid, 'missing_value', missing_value) )
       else
          missing_value=0.0
       endif

       do member = 1, dim_ens

          if (verbose>0 .and. mype==0 .and. i==1) &
               write (*,'(a,4x,a,i6)') 'NEMO-PDAF','--- read member', member

          if (sfields(i)%ndims == 3) then
             call check( nf90_get_var(ncid, varid, tmp_4d, &
                  start=(/istart, jstart, 1, member/), count=(/ni_p, nj_p, nlvls, 1/)) )
          else
             call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
                  start=(/istart, jstart, member/), count=(/ni_p, nj_p, 1/)) )
          end if

          ! Convert field to state vector
          call field2state(tmp_4d, ens(:,member), sfields(i)%off, sfields(i)%ndims, missing_value)

       enddo

       call check( nf90_close(ncid) )

    end do

    do member = 1, dim_ens
       ! Potentially transform fields
       call transform_field_mv(1, ens(:,member))
    end do

    if (verbose>2) then
       do i = 1, n_fields
          write(*,'(a, 1x, a, a10, 1x, a,1x, 2es13.6)') &
               'NEMO-PDAF','Ensemble min and max for ',trim(sfields(i)%variable),' :     ', &
               minval(ens(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim,:)), &
               maxval(ens(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim,:))
       enddo
    end if

  end subroutine read_ens_mv_loop


!===============================================================================  

!!> Read ensemble as state vectors from ensemble file
!!  
  subroutine read_ens(ensfile_fullname, dim_state, dim_ens, ens)

    use netcdf
    
    implicit none
    
    character(len=*), intent(in)    :: ensfile_fullname        !< Name and path of ensemble file
    integer(4),       intent(in)    :: dim_state               !< PE-local state dimension
    integer(4),       intent(in)    :: dim_ens                 !< Ensemble size
    real(8),          intent(inout) :: ens(dim_state, dim_ens) !< Ensemble array
    
    ! Local variables  
    integer(4) :: i                     ! Counter
    integer(4) :: ncid                  ! NC file id
    integer(4) :: dim_state_file        ! state dimension in file
    integer(4) :: dim_ens_file          ! Ensemble size in file
    character(len=400) :: varstr        ! String describing variables in state vector
    character(len=400) :: varstr_file   ! String describing variables in state vector
    integer(4) :: dimstate_dimid, dimens_dimid, ens_varid


! *** Generate string describing the state vector ***
    varstr = ''
    do i = 1, n_fields
       if (i==1) then
          varstr = trim(sfields(i)%variable)
       else
          varstr = trim(varstr)//' '//trim(sfields(i)%variable)
       endif
    end do

! *** Read file
    
    if (verbose>0 .and. mype==0) then
       write(*,'(1x,a,a)') "--- Read ensemble file: ", trim(ensfile_fullname)
    end if
    
    ! Open the file
    call check( nf90_open(ensfile_fullname, nf90_nowrite, ncid) )

    ! Read the string describing the state vector
    call check( nf90_get_att(ncid, NF90_GLOBAL, "state_fields", varstr_file) )

    ! Check consistency of state vector setup
    if (trim(varstr) == trim(varstr_file)) then
    
       ! Get the dimensions
       call check( nf90_inq_dimid(ncid, 'dim_state', dimstate_dimid) )  
       call check(nf90_inquire_dimension(ncid,dimstate_dimid,len=dim_state_file))
  
       call check( nf90_inq_dimid(ncid, 'dim_ens', dimens_dimid) )  
       call check(nf90_inquire_dimension(ncid,dimens_dimid,len=dim_ens_file))

       ! Check consistency of state dimension
       if (dim_state_file == dim_state) then

          ! Check consistency of ensemble size
          if (dim_ens_file >= dim_ens) then

             !  Read ensemble
             call check( nf90_inq_varid(ncid, 'ensemble', ens_varid) )
 
             call check( nf90_get_var(ncid, ens_varid, ens, start=(/1,1/),count=(/dim_state,dim_ens/)) )

             if (dim_ens_file> dim_ens) &
                  write (*,*) 'Notice: Ensemble in file is larger than dim_ens'
          else
             write (*,'(1x,a)') 'ERROR: Ensemble in file is too small'
             write (*,'(1x,a)')  'Stopping program!' 
             stop 10
          end if

       else
          write (*,'(1x,a)') 'ERROR: inconsistent state dimension'
          write (*,'(1x,a)')  'Stopping program!' 
          stop 10
       end if

    else
       write (*,'(1x,a)') 'ERROR: inconsistent variables in state'
       write (*,'(1x,a)')  'Stopping program!' 
       stop 10
    end if

    call check( nf90_close(ncid) )

  end subroutine read_ens


!================================================================================

!> Initialie ensemble array from a list of NEMO output files
!!
  subroutine gen_ens_mv(flate, infilelist, inpath, dim_p, dim_ens, ens)

  implicit none
  
! *** Arguments ***
  real(8),            intent(in)    :: flate        !< inflation
  character(len=*),   intent(in)    :: infilelist   !< Name of file holding dates of input files
  character(len=*),   intent(in)    :: inpath       !< Path to input files
  integer(4),         intent(in)    :: dim_p        !< State dimension
  integer(4),         intent(in)    :: dim_ens      !< Ensemble size
  real(8),            intent(inout) :: ens(:, :)    !< Ensemble array

! *** Local variables ***
  integer(4)        :: i, k, iens, ifile  ! Counters
  integer           :: ios                ! Flag for file reading
  character(len=8)  :: indate             ! Date string
  real(8)           :: ens_mean           ! Ensemble mean
  real(8)           :: invsteps           ! Inverse of ensemble size


! *** Read ensemble from files ***

  if (verbose>0 .and. mype==0) &
       write(*,'(/1x,a)') "*** Generating ensemble from output files ***"

  open (unit=10,file=trim(infilelist),iostat=ios)
  if (ios /= 0) write(*,*) 'Could not open file ',infilelist 
  
  iens=0

  ensloop: do

     read (10,*,iostat=ios) indate
     if (ios/=0) exit ensloop

     do k =1, ntimec
        iens = iens + 1
        if (verbose>0 .and. mype==0) write (*,*) '--- Read ensemble member', iens
        call read_ens_mv(inpath, indate, indate, dim_p, k, ens(:,iens))
     enddo

     if (iens==dim_ens) exit ensloop

  enddo ensloop
  
  close(10)

  ! Check ensemble size
  if (iens<dim_ens) then
     write (*,'(/1x,a)') 'ERROR: Available files less than ensemble size!'
     write (*,'(1x,a)')  'Stopping program!' 
     stop 10
  end if

 
! *** Subtract ensemble mean and inflate ensemble perturbations ***

  invsteps = 1.0/real(dim_ens)


!$OMP PARALLEL DO private(k, ens_mean)
  do k=1,dim_p
     ens_mean = 0.0
     do i=1,dim_ens
        ens_mean = ens_mean + invsteps*ens(k,i)
     end do

     do i=1,dim_ens
        ens(k,i) = flate*(ens(k,i)-ens_mean)
     end do
  end do
!$OMP END PARALLEL DO


end subroutine gen_ens_mv


!=============================================================================== 

!> Read a model field into the state vector of an ensemble array
!!  
  subroutine read_ens_mv(path, date1, date2, dim_p, itime, state)

    use netcdf
    
    implicit none

! *** Arguments ***    
    character(len = *), intent(in)   :: path      !< Path of file
    character(len = *), intent(in)   :: date1     !< First date in file name
    character(len = *), intent(in)   :: date2     !< Second date in file name
    integer(4),         intent(in)   :: dim_p     !< State dimension
    integer(4),         intent(in)   :: itime     !< Time in file to read
    real(8),            intent(inout):: state(:)  !< State vector
    
! *** Local variables ***
    integer(4) :: i, cnt           ! Counters
    integer(4) :: ncid             ! NC file ID
    integer(4) :: varid            ! Variable ID
    character(len=50) :: filename  ! Full file name
    character(len=17) :: dates     ! Combined date string of file

    if (verbose>1) &
         write(*,*) "*** Ensemble: Read model output at time step: ", itime

    ! Initialize state
    state = 0.0

    do i = 1, n_fields

       if (datestype=='_') then
          dates = trim(date1)//'_'//trim(date2)
       else
          dates = trim(date1)//'-'//trim(date2)
       end if
       filename = trim(sfields(i)%file)//trim(dates)//trim(sfields(i)%file_post)//'.nc'

       if (verbose>1) then
          write(*,*) trim(path)//trim(filename)
          write (*,*) i, ' Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       ! Open the file
       call check( nf90_open(trim(path)//trim(filename), nf90_nowrite, ncid) )

       !  Read field
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), varid) )

       if (sfields(i)%ndims == 3) then
          call check( nf90_get_var(ncid, varid, tmp_4d, &
               start=(/istart, jstart, 1, itime/), count=(/ni_p, nj_p, nlvls, 1/)) )
       else
          call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
               start=(/istart, jstart, itime/), count=(/ni_p, nj_p, 1/)) )
       end if

       call check( nf90_close(ncid) )


       ! Convert field to state vector
       call field2state(tmp_4d, state, sfields(i)%off, sfields(i)%ndims, missing_value)

    end do
 
    ! Potentially transform fields
    call transform_field_mv(1, state)

    if (verbose>1) then
       do i = 1, n_fields
          write(*,*) 'Min and max for ',trim(sfields(i)%variable),' :     ',              &
               minval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim)), &
               maxval(state(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim))
       enddo
    end if

  end subroutine read_ens_mv


!================================================================================

!> Initialize array holding model trajectory from a list of NEMO output files
!!
!! Note: This routine is nearly identical to gen_ens_mv, except that the final
!! subtraction of the mean state is not present here.
  subroutine read_trajectory_mv(flate, infilelist, inpath, dim_p, dim_ens, ens)

  implicit none
  
! *** Arguments ***
  real(8),            intent(in)    :: flate        !< inflation
  character(len=*),   intent(in)    :: infilelist   !< Name of file holding dates of input files
  character(len=*),   intent(in)    :: inpath       !< Path to input files
  integer(4),         intent(in)    :: dim_p        !< State dimension
  integer(4),         intent(in)    :: dim_ens      !< Ensemble size
  real(8),            intent(inout) :: ens(:, :)    !< Ensemble array

! *** Local variables ***
  integer(4)        :: i, k, iens, ifile  ! Counters
  integer           :: ios                ! Flag for file reading
  character(len=8)  :: indate             ! Date string
  real(8)           :: ens_mean           ! Ensemble mean
  real(8)           :: invsteps           ! Inverse of ensemble size


! *** Read ensemble from files ***

  if (verbose>0 .and. mype==0) &
       write(*,'(/1x,a)') "*** Read trajectory from output files ***"

  open (unit=10,file=trim(infilelist),iostat=ios)
  if (ios /= 0) write(*,*) 'Could not open file ',infilelist 
  
  iens=0

  ensloop: do

     read (10,*,iostat=ios) indate
     if (ios/=0) exit ensloop

     do k =1, ntimec
        iens = iens + 1
        if (verbose>0 .and. mype==0) write (*,*) '--- Read trajectory state', iens
        call read_ens_mv(inpath, indate, indate, dim_p, k, ens(:,iens))
     enddo

     if (iens==dim_ens) exit ensloop

  enddo ensloop
  
  close(10)

  ! Check ensemble size
  if (iens<dim_ens) then
     write (*,'(/1x,a)') 'ERROR: Available files less than ensemble size!'
     write (*,'(1x,a)')  'Stopping program!' 
     stop 10
  end if

end subroutine read_trajectory_mv

!================================================================================

!> Write an ensemble file holding the state vectors
!!
  subroutine write_state_ens(path, file, dim_state, dim_ens, ens)

    use netcdf

    implicit none 

! *** Arguments *** 
    character(len=*), intent(in):: path          !< Path of file
    character(len=*), intent(in):: file          !< File name
    integer(4),       intent(in):: dim_state     !< state dimension
    integer(4),       intent(in):: dim_ens       !< Ensemble size
    real(8),          intent(in):: ens(:,:)      !< Ensemble array

! *** Local variables *** 
    integer(4) :: i          ! Counter
    integer(4) :: fileid     ! NC file id
    integer(4) :: dimids(2)  ! dimension ids
    integer(4) :: id_ens     ! variable id
    real(8)    :: fillval    ! fill value
    integer(4) :: startv(2),countv(2)  ! Arrays for writing
    character(len=400) :: varstr   ! String describing variables in state vector
    character(len=200) :: filestr  ! String for file name

! *** Generate string describing the state vector ***
    varstr = ''
    do i = 1, n_fields
       if (i==1) then
          varstr = trim(sfields(i)%variable)
       else
          varstr = trim(varstr)//' '//trim(sfields(i)%variable)
       endif
    end do

    if (npes==1) then
       filestr = trim(file)//'.nc'
    else
       filestr = trim(file)//'_'//trim(str(mype))//'.nc'
    end if

! *** Write ensemble of state vectors ***

    ! *** Open file and initialize dimensions and fields *** 
    call check( NF90_CREATE(trim(filestr),NF90_NETCDF4,fileid) )
    call check( NF90_PUT_ATT(fileid,NF90_GLOBAL,'title', &
         'Ensemble matrix for NEMO') )
    call check( nf90_put_att(fileid, NF90_GLOBAL, "state_fields", trim(varstr)) )
  
    ! define dimensions
    call check( NF90_DEF_DIM(fileid,'dim_state',dim_state,dimids(1)) )
    call check( NF90_DEF_DIM(fileid,'dim_ens',dim_ens,dimids(2)) )

    ! define variables
    call check( NF90_DEF_VAR(fileid,'ensemble',NF90_DOUBLE,dimids(1:2),id_ens) )
    fillval = 0.0
    call check( nf90_put_att(fileid, id_ens, "_FillValue", fillval) )
    call check( nf90_put_att(fileid, id_ens, "missing_value", fillval) )
    call check( NF90_def_var_deflate(fileid,id_ens,0,1,1) )

    ! End define mode
    call check( NF90_ENDDEF(fileid) )
 
    do i=1,dim_ens
       startv(1) = 1
       startv(2) = i
       countv(1) = dim_state
       countv(2) = 1
       call check( nf90_put_var(fileid,id_ens,ens(1:dim_state,i), startv,countv) )
    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(fileid) )

  end subroutine write_state_ens


!================================================================================

!> Write ensemble as single files holding model fields
!!
  subroutine write_ens_files(path,file_ens,dim_state,dim_ens,ens)

    use netcdf

    implicit none 
 
! *** Arguments ***
    character(len=*), intent(in) :: path        !< Path of file
    character(len=*), intent(in) :: file_ens    !< Name stub of file
    integer(4),       intent(in) :: dim_state   !< State dimension
    integer(4),       intent(in) :: dim_ens     !< Ensemble size
    real(8),       intent(inout) :: ens(:, :)   !< Ensemble array

! *** Local variables ***
    integer(4)          :: i                ! Counter
    character(len=200)  :: file_ensemble   ! Full name of an ensemble size
    character(len=200)  :: titleEns         ! NC title of file
    real(8)             :: time             ! Time in file


! *** Write ensemble perturbation files

    time=10.0 !TO DO: this is random, time has to be read in pdaf.nml and set here

    titleEns='Ensemble perturbation (ens-mean) for PDAF'
    
    do i=1,dim_ens

       file_ensemble=trim(path)//trim(file_ens)//'_'//trim(str(i))//'.nc'

       call write_field_mv(ens(:, i), dim_state, file_ensemble, titleEns,time, 1, 1)
    enddo

  end subroutine write_ens_files

!================================================================================

!> Write a state vector as model fields into a file
!!
  subroutine write_field_mv(state, dim, filename, title, &
       attime, nsteps, step)
    
    use netcdf
   
    implicit none

! *** Arguments ***
    real(8),          intent(inout) :: state(:)  ! State vector
    integer(4),       intent(in) :: dim          ! State dimension
    character(len=*), intent(in) :: filename     ! File name
    character(len=*), intent(in) :: title        ! File title
    real(8),          intent(in) :: attime       ! Time attribute
    integer(4),       intent(in) :: nsteps       ! Number of time steps stored in file
    integer(4),       intent(in) :: step         ! Time index to write at

! *** Local variables ***
    integer(4) :: ncid
    integer(4) :: dimids_field(4)
    integer(4) :: cnt,i,j,k
    integer(4) :: dimid_time, dimid_lvls, dimid_lat, dimid_lon, dimid_one
    integer(4) :: id_lat, id_lon, id_lev, id_time, id_field
    integer(4) :: dimids(4)
    integer(4) :: startC(2), countC(2)
    integer(4) :: startt(4), countt(4)
    integer(4) :: startz(1), countz(1)
    real(8)    :: fillval
    real(8)    :: timeField(1)

    timeField(1)=attime

    if (step==1) then

! *** Create file ***

       if (verbose>0 .and. mype==0) &
            write (*,'(1x,a,a)') 'Create file: ', trim(filename)

       if (npes==1) then
          call check( NF90_CREATE(trim(filename),NF90_NETCDF4,ncid))
       else
          call check( NF90_CREATE_PAR(trim(filename), NF90_NETCDF4, comm_filter, MPI_INFO_NULL, ncid))
       end if
       call check( NF90_PUT_ATT(ncid, NF90_GLOBAL, 'title', trim(title)))
     
       ! define dimensions for NEMO-input file
       call check( NF90_DEF_DIM(ncid,'t', nsteps, dimid_time))
       call check( NF90_DEF_DIM(ncid, 'z', nlvls, dimid_lvls))
       call check( NF90_DEF_DIM(ncid, 'y', nlats, dimid_lat) )
       call check( NF90_DEF_DIM(ncid, 'x', nlons, dimid_lon) )
       call check( NF90_DEF_DIM(ncid, 'one', 1, dimid_one) )
   
       dimids_field(4)=dimid_time
       dimids_field(3)=dimid_lvls
       dimids_field(2)=dimid_lat
       dimids_field(1)=dimid_lon
       
       ! define variables
       call check( NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, id_time))
       call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
       call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
       call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))
       if (do_deflate) then
          call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
          call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
          call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
       end if

       do i = 1, n_fields
          if (sfields(i)%ndims==3) then
             dimids_field(3)=dimid_lvls
             call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable), NF90_DOUBLE, dimids_field(1:4), id_field) )
          else
             dimids_field(3)=dimid_time
             call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable), NF90_DOUBLE, dimids_field(1:3), id_field) )
          end if
          if (do_deflate) &
               call check( NF90_def_var_deflate(ncid, id_field, 0, 1, 1) )
          call check( nf90_put_att(ncid, id_field, "coordinates", "nav_lat nav_lon") )
          fillval = 1.0e20
          call check( nf90_put_att(ncid, id_field, "_FillValue", fillval) )
          call check( nf90_put_att(ncid, id_field, "missing_value", fillval) )
       end do
       
       ! End define mode
       call check( NF90_ENDDEF(ncid) )

       ! write coordinates
       startz(1)=1
       countz(1)=nlvls

       startC(1)=1
       startC(2)=1
       countC(1)=nlons
       countC(2)=nlats

       if (mype==0) then
          call check( nf90_put_var(ncid,id_lev,depths,startz,countz))
          call check( nf90_put_var(ncid,id_lon,lons,startC,countC))
          call check( nf90_put_var(ncid,id_lat,lats,startC,countC))
       end if

    else
       if (verbose>0) &
            write (*,'(1x,a,a)') 'Open file: ', trim(filename)

       if (npes==1) then
          call check( nf90_open(trim(filename), NF90_WRITE, ncid) )
       else
          call check( nf90_open_par(trim(filename), NF90_WRITE, comm_filter, MPI_INFO_NULL, ncid) )
       end if

    end if


    ! *** Write fields

    call check( nf90_inq_varid(ncid, 'time', id_time) )
!    call check( nf90_VAR_PAR_ACCESS(NCID, id_time, NF90_COLLECTIVE) )
!    call check( nf90_put_vara(ncid, id_time, timeField, start=(/step/), count=(/1/)))
    startt(1) = step
    countt(1) = 1 
    call check( nf90_put_var(ncid, id_time, timeField, startt(1:1), countt(1:1)))

    ! Backwards transformation of state
    call transform_field_mv(2, state)

    do i = 1, n_fields

       ! Convert state vector to field
       tmp_4d = 1.0e20
       call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)

       if (verbose>1) &
            write (*,'(5x,a,a)') '--- write variable: ', trim(sfields(i)%variable)
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable), id_field) )
!       call check( nf90_VAR_PAR_ACCESS(NCID, id_field, NF90_COLLECTIVE) )

       ! Attention with coordinates, in Nemo Restart it is var(time,depth,y,x)
       startt(1) = istart
       countt(1) = ni_p
       startt(2) = jstart
       countt(2) = nj_p
       startt(3) = 1
       countt(3) = nlvls
       startt(4) = step
       countt(4) = 1 
    
       if (sfields(i)%ndims==3) then
          startt(3) = 1
          countt(3) = nlvls
          call check( nf90_put_var(ncid, id_field, tmp_4d, startt, countt))
       else
          startt(3) = step
          countt(3) = 1 

          call check( nf90_put_var(ncid, id_field, tmp_4d, startt(1:3), countt(1:3)))
       end if

    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(ncid) )

  end subroutine write_field_mv

!================================================================================

!> Write a covariance matrix into file
!!
  subroutine write_covar_mv(path, name, dim_p, rank, svdU, meanstate, svals)
    
    use netcdf
   
    implicit none

! *** Arguments ***
    character(len=*), intent(in) :: path         ! File path
    character(len=*), intent(in) :: name         ! File name
    integer(4),       intent(in) :: dim_p        ! State dimension
    integer(4),       intent(in) :: rank         ! Number of EOFs to write
    real(8),          intent(inout) :: svdU(:,:) ! Array of singular vectors
    real(8),          intent(in) :: meanstate(:) ! State vector
    real(8),          intent(in) :: svals(:)     ! Vector of singular values

! *** Local variables ***
    integer(4) :: ncid
    integer(4) :: dimids_field(4)
    integer(4) :: cnt,i,j,k, member
    integer(4) :: dimid_rank, dimid_lvls, dimid_lat, dimid_lon, dimid_one, dimid_state
    integer(4) :: id_lat, id_lon, id_lev, id_time, id_field, id_svals
    integer(4) :: dimids(4)
    integer(4) :: startC(2), countC(2)
    integer(4) :: startt(4), countt(4)
    integer(4) :: startz(1), countz(1)
    !real(8)    :: fillval
    real(4)    :: fillval
    real(8)    :: timeField(1)
    character(len=200) :: filename


! *** Create file ***

!    filename = trim(path)//trim(name)

    if (verbose>0 .and. mype==0) then
       write (*,'(1x,a)') 'Write covariance matrix into file'
       write (*,'(1x,a,a)') 'Create file: ', trim(path)//trim(name)
       if (do_deflate) write (*,'(1x,a)') '--- Apply deflation'
    endif

    if (npes==1) then
       call check( NF90_CREATE(trim(path)//trim(name),NF90_NETCDF4,ncid))
    else
       call check( NF90_CREATE_PAR(trim(path)//trim(name), NF90_NETCDF4, comm_filter, MPI_INFO_NULL, ncid))
    end if
    call check( NF90_PUT_ATT(ncid, NF90_GLOBAL, 'title', 'Covariance matrix'))
    
    ! define dimensions for NEMO-input file
    call check( NF90_DEF_DIM(ncid,'rank', rank, dimid_rank))
    call check( NF90_DEF_DIM(ncid, 'z', nlvls, dimid_lvls))
    call check( NF90_DEF_DIM(ncid, 'y', nlats, dimid_lat) )
    call check( NF90_DEF_DIM(ncid, 'x', nlons, dimid_lon) )
    call check( NF90_DEF_DIM(ncid, 'one', 1, dimid_one) )
    call check( NF90_DEF_DIM(ncid, 'dim_state', dim_p, dimid_state) )
   
    dimids_field(4)=dimid_rank
    dimids_field(3)=dimid_lvls
    dimids_field(2)=dimid_lat
    dimids_field(1)=dimid_lon
       
    ! define variables
    call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
    call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
    call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))

    if (do_deflate) then
       call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
       call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
       call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
    end if

    call check( NF90_DEF_VAR(ncid, 'sigma', NF90_FLOAT, dimids_field(4), id_svals))

    do i = 1, n_fields
       if (sfields(i)%ndims==3) then
          dimids_field(3)=dimid_lvls
          call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable)//'_svd', NF90_FLOAT, dimids_field(1:4), id_field) )
       else
          dimids_field(3)=dimid_rank
          call check( NF90_DEF_VAR(ncid, trim(sfields(i)%variable)//'_svd', NF90_FLOAT, dimids_field(1:3), id_field) )
       end if
       if (do_deflate) then
          call check( NF90_def_var_deflate(ncid, id_field, 0, 1, 1) )
       end if
       call check( nf90_put_att(ncid, id_field, "coordinates", "nav_lat nav_lon") )
       fillval = 1.0e20
       call check( nf90_put_att(ncid, id_field, "_FillValue", fillval) )
       call check( nf90_put_att(ncid, id_field, "missing_value", fillval) )
    end do
       
    ! End define mode
    call check( NF90_ENDDEF(ncid) )

    ! write coordinates
    startz(1)=1
    countz(1)=nlvls

    startC(1)=1
    startC(2)=1
    countC(1)=nlons
    countC(2)=nlats

    if (mype==0) then
       call check( nf90_put_var(ncid,id_lev,depths,startz,countz))
       call check( nf90_put_var(ncid,id_lon,lons,startC,countC))
       call check( nf90_put_var(ncid,id_lat,lats,startC,countC))
    end if

     ! *** Write singular values

    if (mype==0) then
       startz(1)=1
       countz(1)=rank
       call check( nf90_put_var(ncid, id_svals, svals(1:rank), startz, countz))
    end if

    ! *** Write singular vectors as fields

    do member = 1, rank

       ! Backwards transformation of state
       call transform_field_mv(2, svdU(:,member))

    enddo

    do i = 1, n_fields

       if (verbose>1 .and. mype==0) then
          write (*,'(a,i5,a,a,a,i10)') &
               'NEMO-PDAF', i, 'Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       do member = 1, rank

          ! Convert state vector to field
          tmp_4d = 1.0e20
          call state2field(svdU(:,member), tmp_4d, sfields(i)%off, sfields(i)%ndims)

          if (verbose>0 .and. mype==0 .and. i==1) &
               write (*,'(a,4x,a,i6)') 'NEMO-PDAF','--- write rank', member

          call check( nf90_inq_varid(ncid, trim(sfields(i)%variable)//'_svd', id_field) )

          ! Attention with coordinates, in Nemo Restart it is var(time,depth,y,x)
          startt(1) = istart
          countt(1) = ni_p
          startt(2) = jstart
          countt(2) = nj_p
          startt(3) = 1
          countt(3) = nlvls
          startt(4) = member
          countt(4) = 1 
    
          if (sfields(i)%ndims==3) then
             startt(3) = 1
             countt(3) = nlvls
             call check( nf90_put_var(ncid, id_field, tmp_4d, startt, countt))
          else
             startt(3) = member
             countt(3) = 1 

             call check( nf90_put_var(ncid, id_field, tmp_4d, startt(1:3), countt(1:3)))
          end if

       end do
    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(ncid) )

  end subroutine write_covar_mv


!=============================================================================== 

!> Read an covariance matrix from file
!!  
  subroutine read_eof_cov_mv(path, name, dim_p, rank, eofV, svals)

    use netcdf
    use mod_parallel_pdaf, only: abort_parallel
    
    implicit none

! *** Arguments ***    
    character(len=*), intent(in) :: path         ! File path
    character(len=*), intent(in) :: name         ! File name
    integer(4),       intent(in) :: dim_p        ! State dimension
    integer(4),       intent(in) :: rank         ! Number of EOFs to write
    real(8),          intent(inout) :: eofV(:,:) ! Array of singular vectors
    real(8),          intent(inout) :: svals(:)  ! Vector of singular values
    
! *** Local variables ***
    integer(4) :: i, cnt, member   ! Counters
    integer(4) :: ncid             ! NC file ID
    integer(4) :: dimid, varid     ! Dimension and variable ID
    character(len=50) :: filename  ! Full file name
    character(len=17) :: dates     ! Combined date string of file
    integer :: dim_file, rank_file

    if (verbose>0 .and. mype==0) &
         write(*,'(a,4x,a)') 'NEMO-PDAF','--- Read covariance matrix'

    if (.not. allocated(tmp_4d)) allocate(tmp_4d(ni_p, nj_p, nk_p, 1))

    ! Initialize singular vector array
    eofV = 0.0

    filename = trim(path)//trim(name)
    if (verbose>0 .and. mype==0) then
       write(*,'(a,6x,a,1x,a)') 'NEMO-PDAF', '--- File: ', trim(filename)
    end if

    ! Open the file
    call check( nf90_open(trim(filename), nf90_nowrite, ncid) )

    ! Read size of state vector
    call check( NF90_INQ_DIMID(ncid, 'dim_state', dimid) )
    call check( NF90_Inquire_dimension(ncid, dimid, len=dim_file) )
    ! Read rank stored in file
    call check( NF90_INQ_DIMID(ncid, 'rank', dimid) )
    call check( NF90_Inquire_dimension(ncid, dimid, len=rank_file) )

    ! Check consistency of dimensions
    checkdim: IF (rank_file < rank) THEN
      ! *** Rank stored in file is smaller than requested EOF rank ***
      WRITE(*, '(a)') 'Error: Rank stored in file is smaller than requested EOF rank'
      call check( NF90_CLOSE(ncid) )
      call abort_parallel()
    END IF checkdim

    !  Read singular values
    call check( nf90_inq_varid(ncid, 'sigma', varid))
    call check( NF90_GET_VAR(ncid, varid, svals) )

    do i = 1, n_fields

       if (verbose>1 .and. mype==0) then
          write (*,'(a,i5,a,a,a,i10)') &
               'NEMO-PDAF', i, ' Variable: ',trim(sfields(i)%variable), ',  offset', sfields(i)%off
       end if

       !  Read field
       call check( nf90_inq_varid(ncid, trim(sfields(i)%variable)//'_svd', varid) )

       do member = 1, rank

          if (verbose>0 .and. mype==0 .and. i==1) &
               write (*,'(a,6x,a,i6)') 'NEMO-PDAF','--- read EOF rank', member

          if (sfields(i)%ndims == 3) then
             call check( nf90_get_var(ncid, varid, tmp_4d, &
                  start=(/istart, jstart, 1, member/), count=(/ni_p, nj_p, nlvls, 1/)) )
          else
             call check( nf90_get_var(ncid, varid, tmp_4d(:,:,1,1), &
                  start=(/istart, jstart, member/), count=(/ni_p, nj_p, 1/)) )
          end if

          ! Convert field to state vector
          call field2state(tmp_4d, eofV(:,member), sfields(i)%off, sfields(i)%ndims, missing_value)

       enddo

    end do

    call check( nf90_close(ncid) )

    if (verbose>2) then
       do i = 1, n_fields
          write(*,'(a, 1x, a, a10, 1x, a,1x, 2es13.6)') &
               'NEMO-PDAF','Sing. vectors min and max for ',trim(sfields(i)%variable),' :     ', &
               minval(eofV(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim,:)), &
               maxval(eofV(sfields(i)%off+1:sfields(i)%off+sfields(i)%dim,:))
       enddo
    end if

  end subroutine read_eof_cov_mv


!================================================================================

!> Write increment file
!!
  subroutine write_increment(state, dim, filename, id_field)
    
    use netcdf
   
    implicit none

! *** Arguments ***
    real(8),       intent(inout) :: state(:)    !< State vector
    integer(4),       intent(in) :: dim         !< State dimension
    character(len=*), intent(in) :: filename    !< File name
    integer(4),       intent(in) :: id_field    !< Id of field in stte vector

! *** Local variables ***
    integer(4) :: cnt, i, j, k
    integer(4) :: ncid
    integer(4) :: dimids_field(4)
    integer(4) :: dimid_time, dimid_lvls, dimid_lat, dimid_lon
    integer(4) :: id_dateb, id_datef
    integer(4) :: id_lat, id_lon, id_lev, id_time, id_incr
    integer(4) :: dimids(4)
    integer(4) :: startC(2), countC(2)
    integer(4) :: startt(4), countt(4)
    integer(4) :: startz(1), countz(1)
    real(8)    :: fillval
    real(8)    :: timeField(1)
    real(8)    :: bgnTimeInterv(1), finTimeInterv(1), timeInIncr(1)

    ! NOTE: This routine currently only writes a single field from the state vector
    !       which is specified by id_field. this must be a 3D field

    !The time values are read in via namelist at the moment
    ! Fixme
    timeInIncr(1)=incrTime !time for direct initialisation in Nemo (time of restart file which is used for adding to increment file)
    bgnTimeInterv(1)=startEnsTime !Start date of interval on which increment is valid (later for time ramp initialisation of  increment)
    finTimeInterv(1)=endEnsTime !End date of interval on which increment is valid (later for time rap init of increment)! 

    if (verbose>0 .and. mype==0) &
         write (*,'(8x,a)') '--- Write increment file'

    if (npes==1) then
       call check( NF90_CREATE(trim(filename), NF90_NETCDF4, ncid))
    else
       call check( NF90_CREATE_PAR(trim(filename), NF90_NETCDF4, comm_filter, MPI_INFO_NULL, ncid))
    end if

    call check( NF90_PUT_ATT(ncid,  NF90_GLOBAL, 'title', &
         'Increment matrix for NEMO-PDAF coupling'))
     
    ! define dimensions for NEMO-input file
    call check( NF90_DEF_DIM(ncid, 't', 1, dimid_time))
    call check( NF90_DEF_DIM(ncid, 'z', nlvls, dimid_lvls))
    call check( NF90_DEF_DIM(ncid, 'y', nlats, dimid_lat) )
    call check( NF90_DEF_DIM(ncid, 'x', nlons, dimid_lon) )
   
    dimids_field(4)=dimid_time
    dimids_field(3)=dimid_lvls
    dimids_field(2)=dimid_lat
    dimids_field(1)=dimid_lon
       
    ! define variables
    call check( NF90_DEF_VAR(ncid, 'nav_lat', NF90_FLOAT, dimids_field(1:2), id_lat))
    call check( NF90_DEF_VAR(ncid, 'nav_lon', NF90_FLOAT, dimids_field(1:2), id_lon))
    call check( NF90_DEF_VAR(ncid, 'nav_lev', NF90_FLOAT, dimids_field(3), id_lev))
    if (do_deflate) then
       call check( NF90_def_var_deflate(ncid, id_lat, 0, 1, 1) )
       call check( NF90_def_var_deflate(ncid, id_lon, 0, 1, 1) )
       call check( NF90_def_var_deflate(ncid, id_lev, 0, 1, 1) )
    end if

    call check( NF90_DEF_VAR(ncid, 'time', NF90_DOUBLE, id_time))
    call check( NF90_DEF_VAR(ncid, 'z_inc_dateb', NF90_DOUBLE, id_dateb))
    call check( NF90_DEF_VAR(ncid, 'z_inc_datef', NF90_DOUBLE, id_datef)) 

    do i = 1, n_fields

       if (sfields(i)%ndims==3) then
          dimids_field(3)=dimid_lvls
          call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_incr), NF90_DOUBLE, dimids_field(1:4), id_incr) )
       else
          dimids_field(3)=dimid_time
          call check( NF90_DEF_VAR(ncid, trim(sfields(i)%name_incr), NF90_DOUBLE, dimids_field(1:3), id_incr) )
       end if

       if (do_deflate) &
            call check( NF90_def_var_deflate(ncid, id_incr, 0, 1, 1) )

       fillval = 0.0
       call check( nf90_put_att(ncid, id_incr, "long_name", trim(sfields(i)%name_incr)//trim('Increment')) )
       call check( nf90_put_att(ncid, id_incr, "units", trim(sfields(i)%unit)) )
       call check( nf90_put_att(ncid, id_incr, "coordinates", "nav_lat nav_lon") )
       call check( nf90_put_att(ncid, id_incr, "_FillValue", fillval) )
       call check( nf90_put_att(ncid, id_incr, "missing_value", fillval) )
    end do
       
    ! End define mode
    call check( NF90_ENDDEF(ncid) )


    ! *** write coordinates ***
    startz(1)=1
    countz(1)=nlvls

    startC(1)=1
    startC(2)=1
    countC(1)=nlons
    countC(2)=nlats

    if (mype==0) then
       call check( nf90_put_var(ncid, id_lev, depths, startz, countz))
       call check( nf90_put_var(ncid, id_lon, lons, startC, countC))
       call check( nf90_put_var(ncid, id_lat, lats, startC, countC))

       call check( nf90_put_var(ncid, id_time, timeInIncr, start=(/1/), count=(/1/)))
       call check( nf90_put_var(ncid, id_dateb, bgnTimeInterv, start=(/1/), count=(/1/)))
       call check( nf90_put_var(ncid, id_datef, finTimeInterv, start=(/1/), count=(/1/)))
    end if

    ! *** Write fields ***

    ! backwards transformation of increment
    call transform_field_mv(2, state)

    do i = 1, n_fields

       ! Convert state vector to field
       tmp_4d = 0.0
       call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)

       if (verbose>1) &
            write (*,'(5x,a,a)') '--- write variable: ', trim(sfields(i)%name_incr)
       call check( nf90_inq_varid(ncid, trim(sfields(i)%name_incr), id_incr) )

       ! Attention with coordinates, in Nemo Restart it is var(time,depth,y,x)
       startt(1) = istart
       countt(1) = ni_p
       startt(2) = jstart
       countt(2) = nj_p
       startt(3) = 1
       countt(3) = nlvls
       startt(4) = 1
       countt(4) = 1 
           
       if (sfields(i)%ndims==3) then
          countt(3) = nlvls
          call check( nf90_put_var(ncid, id_incr, tmp_4d, startt, countt))
       else
          countt(3) = 1 
          call check( nf90_put_var(ncid, id_incr, tmp_4d, startt(1:3), countt(1:3)))
       end if

    end do

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(ncid) )

  end subroutine write_increment


!================================================================================

!> Overwrite the NEMO restart file
!!
  subroutine update_restart_mv(path, file, state, state_tmp)

    use netcdf
   
    implicit none

! *** Arguments ***
    character(len=*), intent(in) :: path         !< Path for restart file
    character(len=*), intent(in) :: file         !< Restart file name
    real(8),       intent(inout) :: state(:)     !< State vector
    real(8),       intent(inout) :: state_tmp(:) !< tmp state vector (used for storage)

! *** Local variables ***
    integer :: i                      ! Counter
    integer :: ncid                   ! NC file ID
    integer :: lid, uid               ! index range in state vector
    integer :: id_var_n, id_var_b     ! NC variable IDs
    integer  :: startt(4), countt(4)  ! arrays for file writing
    character(len=30) :: rst_file     ! Name of restart file


    ! Attention in run script copy restart file from time of DA to file 'restart_trc_in_befDA.nc' 

    !Write oxy to TRNOXY of restart file (now, time t) -> 
    !restart Nemo with nn_euler=0 (TRBOXY is oxy for t-Delta t)

    ! Store name of restart file
    rst_file = sfields(1)%rst_file

    if (verbose>0) &
         write (*,'(8x,a,1x,a)') '--- Overwrite restart file:',trim(path_restart)//trim(rst_file) 

    ! Open file and retrieve field ids
    if (npes==1) then
       call check( nf90_open(trim(path_restart)//trim(rst_file),NF90_WRITE, ncid))
    else
       call check( nf90_open_par(trim(path_restart)//trim(rst_file),NF90_WRITE,comm_filter, MPI_INFO_NULL, ncid))
    end if

    ! field transformation (if saveIncr=.true. this was already done)
    if (.not.saveIncr) call transform_field_mv(2, state)

    do i = 1, n_fields

       if (trim(sfields(i)%rst_file) /= trim(rst_file)) then
       ! Open other restart file and retrieve field ids
          if (verbose>0 .and. mype==0) &
               write (*,'(8x,a,1x,a)') '--- Open restart file:',trim(path_restart)//trim(sfields(i)%rst_file) 
          if (npes==1) then
             call check( nf90_open(trim(path_restart)//trim(sfields(i)%rst_file),NF90_WRITE, ncid))
          else
             call check( nf90_open_par(trim(path_restart)//trim(sfields(i)%rst_file), &
                  NF90_WRITE, comm_filter, MPI_INFO_NULL, ncid))
          end if
          
          ! Store name of restart file
          rst_file = sfields(i)%rst_file

       end if

       ! Retrieve field IDs
       call check( nf90_inq_varid(ncid, trim(sfields(i)%name_rest_n), id_var_n))
       call check( nf90_inq_varid(ncid, trim(sfields(i)%name_rest_b), id_var_b)) 


       ! backwards transformation state - only if not done by write_increment before

       ! Convert state vector to field
       tmp_4d = 0.0
       call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)

       ! *** write variable for current time ***
       startt(1) = istart
       countt(1) = ni_p
       startt(2) = jstart
       countt(2) = nj_p
       startt(3) = 1
       countt(3) = nlvls
       startt(4) = 1
       countt(4) = 1 
       
       if (sfields(i)%ndims==3) then
          call check( nf90_put_var(ncid, id_var_n, tmp_4d, startt, countt))
       else
          countt(3) = 1
          call check( nf90_put_var(ncid, id_var_n, tmp_4d, startt(1:3), countt(1:3)))
       end if

       ! *** For second (past) time use increment ***

       ! Read field, add increment, and write field
       if (sfields(i)%ndims==3) then
          countt(3) = nlvls
          call check( nf90_get_var(ncid, id_var_b, tmp_4d, startt, countt))
       else
          countt(3) = 1
          call check( nf90_get_var(ncid, id_var_b, tmp_4d, startt(1:3), countt(1:3)))
       end if

       call field2state(tmp_4d, state, sfields(i)%off, sfields(i)%ndims, missing_value)

       lid = sfields(i)%off+1
       uid = sfields(i)%off+sfields(i)%dim
       state(lid : uid) = state(lid : uid) + state_tmp(lid : uid)

       call state2field(state, tmp_4d, sfields(i)%off, sfields(i)%ndims)

       if (sfields(i)%ndims==3) then
          countt(3) = nlvls
          call check( nf90_put_var(ncid, id_var_b, tmp_4d, startt, countt))
       else
          countt(3) = 1
          call check( nf90_put_var(ncid, id_var_b, tmp_4d, startt(1:3), countt(1:3)))
       end if

    end do

    call check( nf90_close(ncid))

  end subroutine update_restart_mv


!================================================================================

!> Check status of NC operation
!!   
  subroutine check(status)

    use netcdf

! *** Aruments ***
    integer, intent ( in) :: status   ! Reading status

    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if

  end subroutine check
 

! ==============================================================================

!> Add a trailing slash to a path string
!!
!! This routine ensures that a string defining a path
!! has a trailing slash.
!!
  subroutine add_slash(path)

    implicit none

! *** Arguments ***
    character(len=100) :: path  !< String holding the path

! *** Local variables ***
    integer :: strlength

! *** Add trailing slash ***
    strlength = len_trim(path)

    if (path(strlength:strlength) /= '/') then
       path = trim(path) // '/'
    end if
    
  end subroutine add_slash


!===============================================================================

!> Convert an integer to a strong of length 4
!!
  character(len=4) function str(k)

    implicit none

    integer, intent(in) :: k   !< number

    write (str, '(i4.4)') k

  end function str


!===============================================================================

!> Check whether a file exists
!!
  function file_exists(filename) result(res)

    implicit none

    character(len=*),intent(in) :: filename   !< File name
    logical                     :: res        !< Status of file

    ! Check if the file exists
    inquire( file=trim(filename), exist=res )

  end function file_exists


!=============================================================================== 
!=== Routines not yet ported
!=============================================================================== 

subroutine gen_ensMeanYearly(flate,varname,infilelist,inpath,dim_state,dim_ens,nyears,transf,trafoConst,ens)

  implicit none
  
  character(len=*),     intent(in)     :: varname
  character(len=80),    intent(in)     :: infilelist
  character(len=160),   intent(in)     :: inpath
  integer(4),           intent(in)     :: dim_state
  integer(4),           intent(in)     :: dim_ens
  integer(4),           intent(in)     :: transf
  real,                 intent(in)     :: trafoConst
  real(8),              intent(inout)  :: ens(dim_state, dim_ens)
  real,                 intent(in)     :: flate
  integer(4),           intent(in)     :: nyears

  real(8)                              :: ens_mean(dim_state)
  real(8)                              :: ens_tmp(dim_state)
  integer(4)                           :: iens,ifile,year
  integer(4)                           :: i,k
  real(8)                              :: invsteps
  character(len=60)                    :: infile
  character(len=220)                   :: infile_full
  integer(4)                           :: ios
  integer(4)                           :: dim_ensOneYear
  real(8), allocatable                 :: ens2(:,:,:) 

  write(*,*) "*** Generating ensemble from model output files with taking mean from every year***"

  dim_ensOneYear=dim_ens/nyears

  allocate(ens2(dim_state,dim_ensOneYear,nyears)) 

  open (unit=10,file=trim(infilelist),iostat=ios)
  if (ios /= 0) write(*,*) 'Could not open file ',infilelist 
  
  iens=0
  yrloop: do year = 1,nyears
    do ifile = 1,int(nfiles/nyears) 
      read (10,*,iostat=ios) infile
      if (ios/=0) exit yrloop
       infile_full=trim(inpath)//trim(infile)
       write(*,*)"file:",infile_full
       do k =1, ntimec
         call read_ergom_ens(varname,infile_full, dim_state, k, transf, trafoConst, ens_tmp)
         iens = iens + 1
         ens2(:,iens,year) = ens_tmp 
       enddo
     enddo
     iens=0
  enddo yrloop
  
  close(10)
 ! write(*,'(/10x,a,i6,i6)') ' dimension of ensemble', dim_ens, iens
 ! write(*,'(/10x,a,i6,i6)') ' dimension of ensemble each year', dim_ensOneYear
  
  invsteps = 1.0/real(dim_ensOneYear)
  iens = 0

  do year=1,nyears
      ens_mean = 0.
     !Parallelisation over spatial component
      !$OMP PARALLEL private(k)
       !$OMP DO
       do k=1,dim_state 
          do i=1,dim_ensOneYear
            ens_mean(k) = ens_mean(k) + invsteps*ens2(k,i,year)
          end do
        enddo
       !$OMP END DO
       !$OMP END PARALLEL
       !write(*,*)'On Master: ens_mean(495450)', ens_mean(495450)!nur zum Testen

        do i=1,dim_ensOneYear
           iens = iens+1
           do k=1,dim_state
  !             write(*,*)'iens',iens
               ens(k,iens) = flate*(ens2(k,i,year)-ens_mean(k))
           end do
        end do
  enddo

  deallocate(ens2)

 end subroutine gen_ensMeanYearly
!================================================================================

!================================================================================
  subroutine gen_ensFlateZ(flate,varname,infilelist,inpath,dim_state,dim_ens, transf, trafoConst, ens)

  implicit none
  
  character(len=3),     intent(in)     :: varname
  character(len=80),    intent(in)     :: infilelist
  character(len=160),   intent(in)     :: inpath
  integer(4),           intent(in)     :: dim_state
  integer(4),           intent(in)     :: dim_ens
  integer(4),           intent(in)     :: transf
  real,                 intent(in)     :: trafoConst
  real(8),              intent(inout)  :: ens(dim_state, dim_ens)
  real,                 intent(in)     :: flate(dim_state)

  real(8)                              :: ens_mean(dim_state)
  real(8)                              :: ens_tmp(dim_state)
  integer(4)                           :: iens,ifile
  integer(4)                           :: i,k
  real(8)                              :: invsteps
  character(len=60)                    :: infile
  character(len=220)                   :: infile_full
  integer(4)                           :: ios

  write(*,*) "*** Generating ensemble from model output files ***"

    open (unit=10,file=trim(infilelist),iostat=ios)
    if (ios /= 0) write(*,*) 'Could not open file ',infilelist 
  
  iens=0
  fileloop: do ifile = 1, nfiles
    read (10,*,iostat=ios) infile
    if (ios/=0) exit fileloop
     infile_full=trim(inpath)//trim(infile)
     write(*,*)"file:",infile_full
     do k =1, ntimec
       call read_ergom_ens(varname,infile_full, dim_state, k, transf, trafoConst, ens_tmp)
       iens = iens + 1
       ens(:,iens) = ens_tmp 
     enddo
  enddo fileloop
  
  close(10)
  write(*,'(/10x,a,i6,i6)') ' dimension of ensemble', dim_ens, iens
  ens_mean = 0.
  invsteps = 1.0/real(dim_ens)

 !Parallelisation over spatial component
  !$OMP PARALLEL private(k)
   !$OMP DO
    do k=1,dim_state 
      do i=1,dim_ens
        ens_mean(k) = ens_mean(k) + invsteps*ens(k,i)
       end do
     enddo

    !$OMP END DO
  !$OMP END PARALLEL

 !write(*,*)'On Master: ens_mean(495450)', ens_mean(495450)!nur zum Testen

  !$OMP PARALLEL private(k)
    !$OMP DO
     do k=1,dim_state
       do i=1,dim_ens
          ens(k,i) = flate(k)*(ens(k,i)-ens_mean(k))
       end do
     end do
     !$OMP END DO
   !$OMP END PARALLEL

  !write(*,*) ens(495450,:)

 end subroutine gen_ensFlateZ
!================================================================================

subroutine read_ens_dims_dim_ens_files(ensfile_path,dim_state,dim_ens)
  
  !! Get the dimension of a varibale from ensemble input file 

    use netcdf
    
    implicit none
    
    character(len = *), intent(in)    :: ensfile_path
    integer(4),           intent(out)   :: dim_state
    integer(4),           intent(out)   :: dim_ens
    
    ! Local variables  
    integer(4)     :: ios
    integer(4)     :: ncid
    logical        :: ens_file_exists
    character(len=220) :: ensfile_fullname
    integer(4) :: dimstate_dimid, dimens_dimid, ens_varid

    dim_ens=0
    ens_file_exists=.true.

    do while(ens_file_exists) 
      ensfile_fullname=trim(ensfile_path)//'ens_'//trim(str(dim_ens+1))//'.nc'
      if (file_exists(ensfile_fullname)) then
         dim_ens=dim_ens+1
      else
         ens_file_exists=.false.
      endif
    enddo

    write(*,*) 'Number of ensemble states: ',dim_ens
    
    ensfile_fullname=trim(ensfile_path)//'ens_'//trim(str(1))//'.nc'

    ! Read ensemble file
    write(*,*) "*** Reading ensemble file to get the dimension: "
    write(*,*) trim(ensfile_fullname)
    
    ! Open the file
    call check( nf90_open(ensfile_fullname, nf90_nowrite, ncid) )
    ! Get the dimensions
    call check( nf90_inq_dimid(ncid, "len", dimstate_dimid) )  
    call check(nf90_inquire_dimension(ncid,dimstate_dimid,len=dim_state))
    ! Close the file. 
    call check( nf90_close(ncid) )
    
   
  end subroutine read_ens_dims_dim_ens_files    
!=============================================================================== 

!================================================================================ 
 subroutine read_ens_dim_ens_files(ensfile_path,dim_state,dim_ens,ens)

 use netcdf

 implicit none

    character(len = *), intent(in)     :: ensfile_path
    integer(4),           intent(in)     :: dim_state               ! PE-local state dimension
    integer(4),           intent(in)     :: dim_ens                 ! 
    real(8),              intent(inout)  :: ens(dim_state, dim_ens)
    
    ! Local variables  
    integer(4) :: i
    integer(4)     :: ios
    integer(4)     :: ncid
    real(8)        :: ens_tmp(dim_state)
    character(len=220) :: ensfile_fullname
    integer(4) :: ens_varid

 
   do i=1,dim_ens
     ensfile_fullname=trim(ensfile_path)//'ens_'//trim(str(i))//'.nc'
 
     write(*,*) "*** Reading ensemble file "
     write(*,*) trim(ensfile_fullname) 

      ! Open the file
      call check( nf90_open(ensfile_fullname, nf90_nowrite, ncid) )
     !  Read ensemble
      call check( nf90_inq_varid(ncid, "ens", ens_varid) )
      call check( nf90_get_var(ncid, ens_varid, ens_tmp, start=(/1/),count=(/dim_state/)) )
      !write(*,*) shape(ens)    
     call check( nf90_close(ncid) )
   
     ens(:,i)=ens_tmp(:)
   enddo

  print*,'ens(495450,:)',ens(495450,:)


 end subroutine read_ens_dim_ens_files
!================================================================================

!================================================================================
   
  subroutine flate_depth(dim_p,nLevFB,nLevFE,flateTOP,flateBOT,flate)

   integer, intent ( in) :: dim_p
   integer, intent(in) :: nLevFB, nLevFE
   real, intent(in) :: flateTOP, flateBOT
   real,intent(out):: flate(dim_p) 
   integer :: x
   integer(4) :: i, j, k


   !make flate factor to deflate ensemble z-dependend
      !ALLOCATE(flate(dim_p)) 
       x = 0
       do k=1,nlvls
         do j = 1,nlats
           do i = 1, nlons
             x=x+1 
             if (k.lt.nLevFB) then 
               flate(x) = flateTOP
             else if ((k.lt.nLevFE).and.(k.gt.nLevFB)) then
               flate(x) = flateTOP-((flateTOP-flateBOT)/(nLevFE-nLevFB))*(k-nLevFB)
             else
               flate(x)=flateBOT
             endif
           enddo
         enddo
       enddo
   
     end subroutine flate_depth


!===============================================================================
! === ROUTINES NO LONGER USED
!===============================================================================

  subroutine read_ens_dims(ensfile_fullname,dim_state,dim_ens)
  
  !! Get the dimension of a variable from ensemble input file 

    use netcdf
    
    implicit none
    
    character(len = *), intent(in)    :: ensfile_fullname
    integer(4),           intent(out)   :: dim_state
    integer(4),           intent(out)   :: dim_ens
    
    ! Local variables  
    integer(4) :: ios
    integer(4) :: ncid
    integer(4) :: dimstate_dimid, dimens_dimid
      
    
    ! Read ensemble file
    write(*,*) "*** Reading ensemble file to get the dimension: "
    write(*,*) trim(ensfile_fullname)
    
    ! Open the file
    call check( nf90_open(ensfile_fullname, nf90_nowrite, ncid) )
    
    ! Get the dimensions
    call check( nf90_inq_dimid(ncid, 'dim_state', dimstate_dimid) )  
    call check(nf90_inquire_dimension(ncid,dimstate_dimid,len=dim_state))
  
    call check( nf90_inq_dimid(ncid, 'dim_ens', dimens_dimid) )  
    call check(nf90_inquire_dimension(ncid,dimens_dimid,len=dim_ens))

    ! Close the file. 
    call check( nf90_close(ncid) )
    
   
  end subroutine read_ens_dims    

!=============================================================================== 

 
  subroutine read_ergom_dims(bio_varname,file_fullname,dim_state_ergom,  &
                        ntimec, nlvls, nlats, nlons, lats, lons, depths, lat1, lon1)
  
  !! Get the dimension of a varibale with using ERGOM netcdf output  

    use netcdf
    
    implicit none
    
    character(len = *),   intent(in)    :: bio_varname
    character(len = *),   intent(in)    :: file_fullname 
    integer(4),           intent(out)   :: dim_state_ergom
    integer(4),           intent(out)   :: ntimec, nlvls, nlats, nlons
    real(8), allocatable, intent(out)   :: lons(:,:), lats(:,:)
    real(8), allocatable, intent(out)   :: depths(:), lat1(:), lon1(:)
    
    ! Local variables  
    integer(4)     :: ios
    integer(4)     :: ncid
    integer(4)     :: i,j  
    integer(4)     :: time_dimid, lvl_dimid, lat_dimid, lon_dimid
    integer(4)     :: depth_varid, bio_varid 
    integer(4)     :: lon_varid, lat_varid
    
    ! Read namelist to get the filename for dimension calculation
    write(*,*) "*** Reading ERGOM output file to get the dimension: "
    write(*,*) trim(file_fullname)
    
    ! Open the file
    call check( nf90_open(trim(file_fullname), nf90_nowrite, ncid) )
    
    ! Get the dimensions
    call check( nf90_inq_dimid(ncid, 'x', lon_dimid) )  
    call check(nf90_inquire_dimension(ncid,lon_dimid,len=nlons))
    
    call check( nf90_inq_dimid(ncid, 'y', lat_dimid) )  
    call check(nf90_inquire_dimension(ncid,lat_dimid,len=nlats))
  
    call check( nf90_inq_dimid(ncid,'deptht', lvl_dimid) )  
    call check(nf90_inquire_dimension(ncid,lvl_dimid,len=nlvls))
  
    call check( nf90_inq_dimid(ncid,'time_counter', time_dimid) )  
    call check(nf90_inquire_dimension(ncid,time_dimid,len=ntimec))   
    
    write(*,*) bio_varname
    ! check if the bio_varname exist
    call check( nf90_inq_varid(ncid, trim(bio_varname), bio_varid) )

    ! calculate the dimension
    dim_state_ergom = nlons*nlats*nlvls

    write(*,*) "       nlons       nlats         nlvls        ntime" 
    write(*,'(i12,i12,i12,i12)') nlons,nlats,nlvls,ntimec
    write(*,'(2x, a, i10)') 'state dimension:', dim_state_ergom
    
    ! read lon, lat and depth 
    allocate(lats(nlons, nlats), lons(nlons, nlats))
    allocate(depths(nlvls))
    
    call check( nf90_inq_varid(ncid, 'nav_lat', lat_varid) )
    call check( nf90_get_var(ncid, lat_varid, lats) )
  
    call check( nf90_inq_varid(ncid, 'nav_lon', lon_varid) )
    call check( nf90_get_var(ncid, lon_varid, lons) )   
    
    call check( nf90_inq_varid(ncid, 'deptht', depth_varid) )
    call check( nf90_get_var(ncid, depth_varid, depths) )

    ! Close the file. 
    call check( nf90_close(ncid) )
   
  end subroutine read_ergom_dims    

! ===================================================================================

  subroutine read_ergom(bio_varname, file_fullname, dim_p, itime, coupling, GaussTransf,trafoConst, state)
  !! read varibale with using ERGOM netcdf output  

    use netcdf
    
    implicit none
    
    character(len = *),   intent(in)   :: bio_varname
    character(len = *), intent(in)   :: file_fullname 
    integer(4),           intent(in)   :: dim_p                ! PE-local state dimension
    integer(4),           intent(in)   :: itime
    character(len = *),   intent(in)   :: coupling
    integer(4),           intent(in)   :: GaussTransf
    real,                 intent(in)   :: trafoConst
    real(8),              intent(inout)  :: state(dim_p)
    
    ! Local variables  
    integer(4)     :: ios
    integer(4)     :: ncid
    integer(4)     :: cnt
    integer(4)     :: time_dimid, lvl_dimid, lat_dimid, lon_dimid
    integer(4)     :: bio_varid


    write(*,*) "*** Reading ERGOM output file at time step: ", itime
    write(*,*) trim(file_fullname)

    print*,'varname ',trim(bio_varname)
    
    ! Open the file
    call check( nf90_open(file_fullname, nf90_nowrite, ncid) )
    
    ! Get the dimensions
    call check( nf90_inq_dimid(ncid, 'x', lon_dimid) )  
    call check(nf90_inquire_dimension(ncid,lon_dimid,len=nlons))
  
    call check( nf90_inq_dimid(ncid, 'y', lat_dimid) )  
    call check(nf90_inquire_dimension(ncid,lat_dimid,len=nlats))
    
    if (coupling=='rest')  then
!       call check( nf90_inq_dimid(ncid,'nav_lev', lvl_dimid) )  
       call check( nf90_inq_dimid(ncid,'deptht', lvl_dimid) )  
       call check(nf90_inquire_dimension(ncid,lvl_dimid,len=nlvls))
    else
       call check( nf90_inq_dimid(ncid,'deptht', lvl_dimid) )  
       call check(nf90_inquire_dimension(ncid,lvl_dimid,len=nlvls))
    endif

    call check( nf90_inq_dimid(ncid,'time_counter', time_dimid) )  
    call check(nf90_inquire_dimension(ncid,time_dimid,len=ntimec))
    
    call check( nf90_inq_varid(ncid, trim(bio_varname), bio_varid) )
    if (coupling=='incr') then
       call check( nf90_get_att(ncid, bio_varid, 'missing_value', missing_value) )
       write(*,*) missing_value
    else
       missing_value=0.0
    endif

    ! Read bio-parameters
    call check( nf90_inq_varid(ncid, trim(bio_varname), bio_varid) )
    
!     call check( nf90_get_var(ncid, bio_varid, bio_values, &
!          start=(/1,1,1,itime/), count=(/nlons,nlats,nlvls,itime/)) )

    call check( nf90_get_var(ncid, bio_varid, tmp_4d, &
         start=(/1,1,1,itime/), count=(/nlons,nlats,nlvls,1/)) )
    !write(*,*) shape(bio_values)   

    call check( nf90_close(ncid) )

write (*,*) 'read_ergom with offset', sfields(id%temp)%off
    ! Convert field to state vector
    call field2state(tmp_4d, state, sfields(id%temp)%off, sfields(id%temp)%ndims, missing_value)
   
!    write(*,*) '5953186->',state(5953186)
   
   write(*,'(1x,a,a10,1x, a,5x, 2f12.6)') &
        'Min and max for ', trim(bio_varname), ':', minval(state), maxval(state)
   
   select case (GaussTransf)
     case(0)
        write(*,*) 'No Transformation of bio variable'
     case(1)
        write(*,*) 'use log basis 10 transformation of bio variable'
        state = log10(state+trafoConst)
     case(2)
        write(*,*) 'use ln transformation of bio variable'
        state = log(state+trafoConst)
     case(3)
        write(*,*) 'no transformation- box cox still needs to be implemented'
     case DEFAULT
        write(*,*) 'No Transformation of bio variable'
   end select

  end subroutine read_ergom

!================================================================================

!================================================================================
 subroutine read_ens_dims_gen_ens(infilelist)

    use netcdf
    
    implicit none
    
    character(len = *),   intent(in)    :: infilelist
    
    ! Local variables  
    integer(4)     :: ios
    integer(4)     :: ncid
    character(80)      :: infile
      
    
    ! Read ensemble file
    write(*,*) "*** Reading ensemble file list to get the dimension: "

   ! Get number of restart files
    nfiles=0
    open (unit=10,file=trim(infilelist),iostat=ios)
    if (ios /= 0) write(*,*) 'Could not open file ',infilelist 

    readloop: do
       read (10,*,iostat=ios) infile
       if (ios/=0) exit readloop
       nfiles=nfiles+1
    enddo readloop
    close(10)
    write(*,'(5x,a,a,a,i6)') ' Number of files listed in ', trim(infilelist),' =', nfiles

   
  end subroutine read_ens_dims_gen_ens    
!================================================================================

!================================================================================
  subroutine gen_ens(flate,varname,infilelist,inpath,dim_state,dim_ens, transf, trafoConst, ens)

  implicit none
  
  character(len=*),     intent(in)     :: varname
  character(len=80),    intent(in)     :: infilelist
  character(len=160),   intent(in)     :: inpath
  integer(4),           intent(in)     :: dim_state
  integer(4),           intent(in)     :: dim_ens
  integer(4),           intent(in)     :: transf
  real,                 intent(in)     :: trafoConst
  real(8),              intent(inout)  :: ens(dim_state, dim_ens)
  real,                 intent(in)     :: flate

  real(8)                              :: ens_mean
  integer(4)                           :: iens,ifile
  integer(4)                           :: i,k
  real(8)                              :: invsteps
  character(len=60)                    :: infile
  character(len=220)                   :: infile_full
  integer(4)                           :: ios

  write(*,'(/1x,a)') "*** Generating ensemble from ERGOM output files ***"

    open (unit=10,file=trim(infilelist),iostat=ios)
    if (ios /= 0) write(*,*) 'Could not open file ',infilelist 
  
  iens=0
  fileloop: do ifile = 1, nfiles
    read (10,*,iostat=ios) infile
    if (ios/=0) exit fileloop
     infile_full=trim(inpath)//trim(infile)
     write(*,*)"ERGOM file:",infile_full
     do k =1, ntimec
       iens = iens + 1
       call read_ergom_ens(varname,infile_full, dim_state, k, transf, trafoConst, ens(:,iens))
     enddo
  enddo fileloop
  
  close(10)

  write(*,'(/10x,a,i6,i6)') ' dimension of ensemble', dim_ens, iens

  invsteps = 1.0/real(dim_ens)

 !Parallelisation over spatial component

!$OMP PARALLEL DO private(k, ens_mean)
  do k=1,dim_state
     ens_mean = 0.0
     do i=1,dim_ens
        ens_mean = ens_mean + invsteps*ens(k,i)
     end do

     do i=1,dim_ens
        ens(k,i) = flate*(ens(k,i)-ens_mean)
     end do
  end do
!$OMP END PARALLEL DO

 end subroutine gen_ens
!================================================================================

!=============================================================================== 
  
  subroutine read_ergom_ens(name_biovar,file_fullname, dim_p, itime, transf, trafoConst, ens)

    use netcdf
    
    implicit none
    
    character(len = *),   intent(in)   :: name_biovar
    character(len = 220), intent(in)   :: file_fullname 
    integer(4),           intent(in)   :: dim_p                ! PE-local state dimension
    integer(4),           intent(in)   :: itime
    integer(4),           intent(in)   :: transf
    real,                 intent(in)   :: trafoConst
    real(8),              intent(inout):: ens(dim_p)
    
    ! Local variables  
    integer(4)     :: ios
    integer(4)     :: ncid
    integer(4)     :: xnan
    integer(4)     :: cnt
    integer(4)     :: bio_varid


    write(*,*) "*** Reading output file at time step: ", itime
    write(*,*) 'File: ', trim(file_fullname)
    
    ! Open the file
    call check( nf90_open(file_fullname, nf90_nowrite, ncid) )

    !  Read bio-parameters
!    allocate( bio_values(nlons, nlats, nlvls,ntimec))
    call check( nf90_inq_varid(ncid, trim(name_biovar), bio_varid) )

    call check( nf90_get_att(ncid, bio_varid, 'missing_value', missing_value) )
!    write(*,*) missing_value

    call check( nf90_get_var(ncid, bio_varid, tmp_4d, start=(/1,1,1,1/),count=(/nlons,nlats,nlvls,ntimec/)) )

    call check( nf90_close(ncid) )


    ! Convert field to state vector

    ens = 0.0
    call field2state(tmp_4d, ens, sfields(id%temp)%off, sfields(id%temp)%ndims, missing_value)
 
    write(*,*) 'Min and max for ',trim(name_biovar),' :     ',              &
         minval(ens), maxval(ens)

    select case (transf)
    case(0)
       write(*,*) 'No Transformation of bio variable'
    case(1)
       write(*,*) 'use log basis 10 transformation of bio variable'
       ens = log10(ens+trafoConst)
    case(2)
       write(*,*) 'use ln transformation of bio variable'
       ens = log(ens+trafoConst)
    case(3)
       write(*,*) 'no transformation- box cox still needs to be implemented'
    case DEFAULT
       write(*,*) 'No Transformation of bio variable'
    end select

  end subroutine read_ergom_ens
!================================================================================

!================================================================================

  subroutine write_nc_field(dummyvector,dim_dummy,filename,fieldname,title,attime)

    use netcdf

    implicit none

    integer(4),         intent(in):: dim_dummy
    real,                intent(in):: dummyvector(dim_dummy)
    character(len=100), intent(in):: filename
    character(len=200), intent(in):: title
    character(len=9),  intent(in):: fieldname
    real,               intent(in):: attime

    integer(4)                    :: fileidField
    integer(4)                    :: dimids_field(4)
    integer(4)                    :: i,j,k, cnt
    integer(4)                    :: idtime,idlvls,idlat,idlon
    integer(4)                    :: id_lat,id_lon,id_lev,id_time,id_field
    integer(4)                    :: startv(4),countv(4),dimids(4),startz(1),countz(1)
    integer(4)                    :: startC(2),countC(2),startt(4),countt(4)
    real                          :: fillval
    integer(4),dimension(10)      :: stat
    real                          :: timeField(1)

    timeField(1)=attime

    ! Convert state vector to field
    tmp_4d = 0.0
    call state2field(dummyvector, tmp_4d, sfields(id%temp)%off, sfields(id%temp)%ndims)

       call check( NF90_CREATE(trim(filename),NF90_NETCDF4,fileidField))
       call check( NF90_PUT_ATT(fileidField,NF90_GLOBAL,'title', &
            title))
     
       ! define dimensions for NEMO-input file
       call check( NF90_DEF_DIM(fileidField,'t',1,idtime))!instead of 1, in normal restart file unlimited
       call check( NF90_DEF_DIM(fileidField,'z',nlvls,idlvls))
       call check( NF90_DEF_DIM(fileidField,'y',nlats,idlat) )
       call check( NF90_DEF_DIM(fileidField,'x',nlons,idlon) )
   
       dimids_field(4)=idtime
       dimids_field(3)=idlvls
       dimids_field(2)=idlat
       dimids_field(1)=idlon
       
       ! define variables
       call check( NF90_DEF_VAR(fileidField,'nav_lat',NF90_FLOAT,dimids_field(1:2),id_lat))!TO DO:How to get float in nc from fortran?  
       call check( NF90_DEF_VAR(fileidField,'nav_lon',NF90_FLOAT,dimids_field(1:2),id_lon)) !TO DO: float variables
       call check( NF90_DEF_VAR(fileidField,'nav_lev',NF90_FLOAT,dimids_field(3),id_lev)) !TO DO: float variables 
       call check( NF90_DEF_VAR(fileidField,'time',NF90_DOUBLE,id_time))
       call check( NF90_DEF_VAR(fileidField,fieldname,NF90_DOUBLE,dimids_field(1:4),id_field) )
       fillval = 0.0
       call check( NF90_def_var_deflate(fileidField,id_field,0,1,1) )
      ! call check( nf90_put_att(fileidField, id_incr, "long_name", "Difference of Time snap shots of model to mean") )
      ! call check( nf90_put_att(fileidField, id_incr, "units", "mmol/m3") )
       call check( nf90_put_att(fileidField, id_field, "coordinates", "nav_lat nav_lon") )
       call check( nf90_put_att(fileidField, id_field, "_FillValue", fillval) )
       call check( nf90_put_att(fileidField, id_field, "missing_value", fillval) )
       
       ! check status flag
       do j=1,5
         if (stat(j).ne.NF90_NOERR) write(*,*) &
           'NetCDF error in dimension definitions, no.',j, &
           ' file ',fileidField
       end do

       ! End define mode
         call check( NF90_ENDDEF(fileidField) )

         !Attention with coordinates, in Nemo Restart it is var(time,depth,y,x)
           startt(1) = 1
           countt(1) = nlons 
           startt(2) = 1
           countt(2) = nlats
           startt(3) = 1
           countt(3) = nlvls
           startt(4) = 1
           countt(4) = 1 
           
           call check( nf90_put_var(fileidField,id_field,tmp_4d,startt,countt))

           startz(1)=1
           countz(1)=nlvls

           startC(1)=1
           startC(2)=1
           countC(1)=nlons
           countC(2)=nlats
           
           call check( nf90_put_var(fileidField,id_lev,depths,startz,countz))
           call check( nf90_put_var(fileidField,id_lon,lons,startC,countC))
           call check( nf90_put_var(fileidField,id_lat,lats,startC,countC))

           call check( nf90_put_var(fileidField,id_time,timeField,start=(/1/),count=(/1/)))

         ! *** close file with state sequence ***
         call check( NF90_CLOSE(fileidField) )

 end subroutine write_nc_field
!================================================================================

!================================================================================

  subroutine write_field(state, dim, filename, fieldname, title, &
       attime, step, offset)
    
    use netcdf
   
    implicit none

    ! *** Arguments ***
    integer(4),       intent(in) :: dim
    real,             intent(in) :: state(:)
    character(len=*), intent(in) :: filename
    character(len=*), intent(in) :: title
    character(len=*), intent(in) :: fieldname
    real,             intent(in) :: attime
    integer(4),       intent(in) :: step        ! time index to write at
    integer(4),       intent(in) :: offset      ! offset in state vector

    ! *** Local variables ***
    integer  :: ncid
    integer  :: dimids_field(4)
    integer  :: cnt,i,j,k
    integer  :: idtime,idlvls,idlat,idlon
    integer  :: id_lat,id_lon,id_lev,id_time,id_field
    integer  :: dimids(4)
    integer  :: startC(2),countC(2)
    integer  :: startt(4),countt(4)
    integer  :: startz(1),countz(1)
    real     :: fillval
    real     :: timeField(1)

    timeField(1)=attime

    ! Convert state vector to field
    tmp_4d = 0.0
    call state2field(state, tmp_4d, sfields(id%temp)%off, sfields(id%temp)%ndims)


    if (step==1) then
       call check( NF90_CREATE(trim(filename),NF90_NETCDF4,ncid))
       call check( NF90_PUT_ATT(ncid,NF90_GLOBAL,'title', trim(title)))
     
       ! define dimensions for NEMO-input file
       call check( NF90_DEF_DIM(ncid,'t',NF90_UNLIMITED,idtime))
       call check( NF90_DEF_DIM(ncid,'z',nlvls,idlvls))
       call check( NF90_DEF_DIM(ncid,'y',nlats,idlat) )
       call check( NF90_DEF_DIM(ncid,'x',nlons,idlon) )
   
       dimids_field(4)=idtime
       dimids_field(3)=idlvls
       dimids_field(2)=idlat
       dimids_field(1)=idlon
       
       ! define variables
       call check( NF90_DEF_VAR(ncid,'nav_lat',NF90_FLOAT,dimids_field(1:2),id_lat))
       call check( NF90_def_var_deflate(ncid,id_lat,0,1,1) )
       call check( NF90_DEF_VAR(ncid,'nav_lon',NF90_FLOAT,dimids_field(1:2),id_lon))
       call check( NF90_def_var_deflate(ncid,id_lon,0,1,1) )
       call check( NF90_DEF_VAR(ncid,'nav_lev',NF90_FLOAT,dimids_field(3),id_lev))
       call check( NF90_def_var_deflate(ncid,id_lev,0,1,1) )
       call check( NF90_DEF_VAR(ncid,'time',NF90_DOUBLE,id_time))
       call check( NF90_DEF_VAR(ncid,trim(fieldname),NF90_DOUBLE,dimids_field(1:4),id_field) )
       call check( NF90_def_var_deflate(ncid,id_field,0,1,1) )
       call check( nf90_put_att(ncid, id_field, "coordinates", "nav_lat nav_lon") )
       fillval = 0.0
       call check( nf90_put_att(ncid, id_field, "_FillValue", fillval) )
       call check( nf90_put_att(ncid, id_field, "missing_value", fillval) )
       
       ! End define mode
       call check( NF90_ENDDEF(ncid) )

       ! write coordinates
       startz(1)=1
       countz(1)=nlvls

       startC(1)=1
       startC(2)=1
       countC(1)=nlons
       countC(2)=nlats
           
       call check( nf90_put_var(ncid,id_lev,depths,startz,countz))
       call check( nf90_put_var(ncid,id_lon,lons,startC,countC))
       call check( nf90_put_var(ncid,id_lat,lats,startC,countC))

    else

       call check( nf90_open(trim(filename), NF90_WRITE, ncid) )
       call check( nf90_inq_varid(ncid, trim(fieldname), id_field) )
       call check( nf90_inq_varid(ncid, 'time', id_time) )

    end if

    ! Attention with coordinates, in Nemo Restart it is var(time,depth,y,x)
    startt(1) = 1
    countt(1) = nlons 
    startt(2) = 1
    countt(2) = nlats
    startt(3) = 1
    countt(3) = nlvls
    startt(4) = step
    countt(4) = 1 
           
    call check( nf90_put_var(ncid, id_field, tmp_4d, startt, countt))


    call check( nf90_put_var(ncid, id_time, timeField, start=(/step/), count=(/1/)))

    ! *** close file with state sequence ***
    call check( NF90_CLOSE(ncid) )

  end subroutine write_field
!================================================================================



end module mod_io_pdaf

! asking if there are the scale_factor and add_offset attributes
! status = nf90_get_att(ncid,varid,"scale_factor",sf)
! IF (status == -43) sf=1.0
! status = nf90_get_att(ncid,varid,"add_offset",ofs)
! IF (status == -43) ofs = 0.0
! ff=sf*var_data+ofs
! call check(nf90_close(ncid))
! PRINT*,'Reading step ',ist
! PRINT*,''
