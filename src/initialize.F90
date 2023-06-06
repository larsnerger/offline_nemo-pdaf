!$Id$
!> Initialize model dimensions
!!
!! Routine to perform initialization of the model information for
!! PDAF. Here, the global size of the model domain, the global size
!! of the model state vector and the sizes for decomposition of the 
!! state vector need to be initialized.
!! Generally, this could also be joined with the routine init_pdaf().
!!
!! !REVISION HISTORY:
!! 2022-02 - Lars Nerger - Initial code
!! Later revisions - see svn log
!!
subroutine initialize()

  use mpi
  use netcdf
  use mod_parallel_pdaf, &
       only: mype=>mype_world, npes=>npes_world, MPIerr, &
       abort_parallel
  use mod_nemo_pdaf, &
       only: path_dims, file_dims, &
       jpiglo, jpjglo, jpk, glamt, gphit, gdept_1d, &
       tmp_4d, dim_2d, lat1, lon1, nwet, nwet3d, &
       wet_pts, idx_wet_2d, idx_nwet, nlev_wet_2d, use_wet_state, &
       nldi, nldj, nlei, nlej, ni_p, nj_p, nk_p, &
       istart, jstart, dim_2d_p, lat1_p, lon1_p, sdim2d, sdim3d
  use mod_io_pdaf, &
       only: check, add_slash
  use mod_memcount_pdaf, &
       only: memcount

  implicit none

! *** Arguments ***
  ! none

! *** local variables *** 
  integer :: screen=1                           ! Verbosity flag
  integer :: i, j, k                            ! Counters
  integer :: cnt, cnt_all, cnt_layers           ! Counters
  integer :: ncid                               ! nc file ID
  integer :: lon_dimid, lat_dimid, lvl_dimid    ! Dimension IDs
  integer :: id_gphit, id_glamt, id_navlev, id_temp  ! variable IDs
  real(8) :: missing_value                           ! missing value
  logical :: have_pdafnml                            ! Flag whether namelist file is present
  character(len=20) :: varname                  ! Name of variable ot read to determine mask
  character(len=1)  :: decomp='y'               ! Direction of decomposition (x, y)
  logical :: read_decomp=.false.                ! Whether to read domain decomposition from file
  character(len=50) :: file_decomp='decomp.txt' ! Name of decomposition file
  integer :: n_domains_lon, n_domains_lat
  integer, allocatable :: dims_lat(:), dims_lon(:)
  integer, allocatable :: nldi_all(:), nldj_all(:) ! first inner index in i/j direction for all PEs
  integer, allocatable :: nlei_all(:), nlej_all(:) ! last inner index in i/j direction for all PEs
  integer :: nx, ny, nz                         ! Size of 3D grid

  namelist /nemo_nml/ screen, path_dims, file_dims, varname, use_wet_state, &
       read_decomp, file_decomp


! **********************
! *** Initialization ***
! **********************

  varname = 'votemper'
  path_dims = './'
  file_dims = 'nemo_output.nc'

! *** Read namelist file for PDAF if the list is there ***

  inquire( FILE='pdaf.nml', EXIST=have_pdafnml ) 
  if (have_pdafnml) then
     open (500,file='pdaf.nml')
     read (500,NML=nemo_nml)
     close (500)
  end if
  call add_slash(path_dims)


! *************************************
! *** Read dimensions of model grid ***
! *************************************

  if (mype==0) then
     write (*,'(a,2x,a)') 'NEMO-PDAF', '*** Reading NEMO dimensions ***'
     write (*,'(a,2x,a,a)') 'NEMO-PDAF', 'File: ',trim(path_dims)//trim(file_dims)
     write (*,'(a,2x,a,a)') 'NEMO-PDAF', 'Variable used for mask: ', trim(varname)
  end if

  ! Open the file
  call check( nf90_open(trim(path_dims)//trim(file_dims), nf90_nowrite, ncid) )
    
  ! Get the dimensions
  call check( nf90_inq_dimid(ncid, 'x', lon_dimid) )  
  call check(nf90_inquire_dimension(ncid,lon_dimid,len=jpiglo))
    
  call check( nf90_inq_dimid(ncid, 'y', lat_dimid) )  
  call check(nf90_inquire_dimension(ncid,lat_dimid,len=jpjglo))
  
  call check( nf90_inq_dimid(ncid,'deptht', lvl_dimid) )  
  call check(nf90_inquire_dimension(ncid,lvl_dimid,len=jpk))

  nx = jpiglo
  ny = jpjglo
  nz = jpk


! ***********************************
! *** Define domain-decomposition ***
! ***********************************

  allocate(nldi_all(0:npes-1))
  allocate(nldj_all(0:npes-1))
  allocate(nlei_all(0:npes-1))
  allocate(nlej_all(0:npes-1))


  if (npes==1 .or. .not.read_decomp) then
     ! Use pre-defined decomposition

     decomp='y'
     if (npes==1) then
        ! No decomposition
        nldi_all(0) = 1
        nldj_all(0) = 1
        nlei_all(0) = nx
        nlej_all(0) = ny
     elseif (npes==2) then
        if (decomp=='y') then
           ! Decomposition in y
           nldi_all(:) = 1
           nlei_all(:) = nx
           nldj_all(0) = 1
           nldj_all(1) = 500
           nlej_all(0) = 500
           nlej_all(1) = ny
        else
           ! Decomposition in x
           nldj_all(:) = 1
           nlej_all(:) = ny
           nldi_all(0) = 1
           nldi_all(1) = 480
           nlei_all(0) = 480
           nlei_all(1) = nx
        end if
     elseif (npes==4) then
        if (decomp=='y') then
           ! Decomposition in y
           nldi_all(:) = 1
           nlei_all(:) = nx
           nldj_all(0) = 1
           nldj_all(1) = 390
           nldj_all(2) = 500
           nldj_all(3) = 610

           nlej_all(0) = 390
           nlej_all(1) = 500
           nlej_all(2) = 610
           nlej_all(3) = ny
        else
           ! Decomposition in x
           nldj_all(:) = 1
           nlej_all(:) = ny
           nldi_all(0) = 1
           nldi_all(1) = 240
           nldi_all(2) = 600
           nldi_all(3) = 910
           nlei_all(0) = 240
           nlei_all(1) = 600
           nlei_all(2) = 910
           nlei_all(3) = nx
        end if
     elseif (npes==6) then
        ! Decomposition in y
        nldi_all(:) = 1
        nlei_all(:) = nx
        nldj_all(0) = 1
        nldj_all(1) = 350
        nldj_all(2) = 420
        nldj_all(3) = 500
        nldj_all(4) = 570
        nldj_all(5) = 670

        nlej_all(0) = 350
        nlej_all(1) = 420
        nlej_all(2) = 500
        nlej_all(3) = 570
        nlej_all(4) = 670
        nlej_all(5) = ny
     elseif (npes==12) then
        if (decomp=='y') then

           ! Decomposition in y
           nldi_all(:) = 1
           nlei_all(:) = nx
           nldj_all(0) = 1
           nldj_all(1) = 250
           nldj_all(2) = 350
           nldj_all(3) = 390
           nldj_all(4) = 425
           nldj_all(5) = 460
           nldj_all(6) = 500
           nldj_all(7) = 540
           nldj_all(8) = 570
           nldj_all(9) = 610
           nldj_all(10) = 670
           nldj_all(11) = 790

           nlej_all(0) = 250
           nlej_all(1) = 350
           nlej_all(2) = 390
           nlej_all(3) = 425
           nlej_all(4) = 460
           nlej_all(5) = 500
           nlej_all(6) = 540
           nlej_all(7) = 570
           nlej_all(8) = 610
           nlej_all(9) = 670
           nlej_all(10) = 790
           nlej_all(11) = ny
        elseif (decomp=='b') then
           ! 2-dimensional decomposition
           nldi_all(0) = 1
           nldi_all(1) = 50
           nldi_all(2) = 50
           nldi_all(3) = 1
           nldi_all(4) = 1
           nldi_all(5) = 1
           nldi_all(6) = 460
           nldi_all(7) = 630
           nldi_all(8) = 630
           nldi_all(9) = 700
           nldi_all(10) = 730
           nldi_all(11) = 730

           nlei_all(0) = 360
           nlei_all(1) = 460
           nlei_all(2) = 460
           nlei_all(3) = 460
           nlei_all(4) = 460
           nlei_all(5) = 460
           nlei_all(6) = 630
           nlei_all(7) = 1050
           nlei_all(8) = 1050
           nlei_all(9) = 1238
           nlei_all(10) = 1238
           nlei_all(11) = 1070

           nldj_all(0) = 1
           nldj_all(1) = 230
           nldj_all(2) = 350
           nldj_all(3) = 430
           nldj_all(4) = 500
           nldj_all(5) = 570

           nldj_all(6) = 300
           nldj_all(7) = 300
           nldj_all(8) = 450
           nldj_all(9) = 575
           nldj_all(10) = 680
           nldj_all(11) = 815

           nlej_all(0) = 230
           nlej_all(1) = 350
           nlej_all(2) = 430
           nlej_all(3) = 500
           nlej_all(4) = 570
           nlej_all(5) = 655

           nlej_all(6) = 700
           nlej_all(7) = 450
           nlej_all(8) = 575
           nlej_all(9) = 680
           nlej_all(10) = 815
           nlej_all(11) = ny
        end if
     end if
  else

     ! Read decomposition file

     allocate(dims_lon(0:npes-1))
     allocate(dims_lat(0:npes-1))

     if (mype==0) write(*,'(/a,2x,a,a)') &
          'NEMO-PDAF','*** Read domain decomposition: ', trim(file_decomp)

     open(11,FILE=trim(file_decomp))
     read(11,*) n_domains_lon, n_domains_lat
     if (mype==0) write (*,'(a,3x,a,2i5)') 'NEMO-PDAF', 'Number of tasks lon/lat: ', n_domains_lon, n_domains_lat

     do k = 0, n_domains_lon*n_domains_lat-1
        read(11, *)  nldj_all(k), nlej_all(k), nldi_all(k), nlei_all(k), dims_lon(k), dims_lat(k)
     end do
     close(11)

     if (n_domains_lon*n_domains_lat /= npes) then
        write (*,*) 'ERROR: number of domains in file inconsistent with number of processes'
        call abort_parallel()
     end if

     deallocate(dims_lon, dims_lat)

  end if

  ! Set sizes for this MPI task
  nldj = nldj_all(mype)
  nldi = nldi_all(mype)
  nlei = nlei_all(mype)
  nlej = nlej_all(mype)

  ! Set local domain size
  ni_p = nlei - nldi + 1
  nj_p = nlej - nldj + 1
  nk_p = jpk

  ! Start indices for sub-domain without halo
  istart = nldi
  jstart = nldj

  ! Screen output
  if (mype==0 .and. npes>1 .and. screen>0) then
     write (*,'(/a,3x,a)') 'NEMO-PDAF','Grid decomposition:' 
     write (*,'(a, 8x,a,2x,a,a,2x,a,a,1x,a,6(1x,a))') &
          'NEMO-PDAF','rank ', 'istart', '  iend', 'jstart', '  jend', '  idim', '  jdim'
     do i = 0, npes-1
        write (*,'(a,2x, a,i6,1x,2i7,2i7,2i7)') 'NEMO-PDAF', 'RANK',i, nldj_all(i), nlej_all(i), &
             nldi_all(i), nlei_all(i), nlei_all(i)-nldi_all(i)+1, nlej_all(i)-nldj_all(i)+1
     end do
  end if


! *************************************************
! *** Read coordinates and field to define mask ***
! *************************************************

  call check( nf90_inq_varid(ncid, 'nav_lat', id_gphit) )
  call check( nf90_inq_varid(ncid, 'nav_lon', id_glamt) )
  call check( nf90_inq_varid(ncid, 'deptht', id_navlev) )
  call check( nf90_inq_varid(ncid, trim(varname), id_temp) )

  call check( nf90_get_att(ncid, id_temp, 'missing_value', missing_value) )

  allocate(glamt(nx, ny), gphit(nx, ny))
  allocate(gdept_1d(nz))
  allocate(tmp_4d(ni_p, nj_p, nz, 1))
  call memcount(1, 'r', 2*nx*ny + nz + ni_p*nj_p*nz)

  call check( nf90_get_var(ncid, id_glamt, glamt(:,:), (/1, 1/), (/nx, ny/) ) )
  call check( nf90_get_var(ncid, id_gphit, gphit(:,:), (/1, 1/), (/nx, ny/) ) )
  call check( nf90_get_var(ncid, id_navlev, gdept_1d(:), (/1/), (/nz/) ) )

  call check( nf90_get_var(ncid, id_temp, tmp_4d(:,:,:,:), (/nldi, nldj, 1, 1/), &
       (/ni_p, nj_p, nz, 1/) ) )
  
  ! Close the file. 
  call check( nf90_close(ncid) )


! *** Initialize coordinate vectors ***

  allocate(lat1(ny), lon1(nx))
  allocate(lat1_p(nj_p), lon1_p(ni_p))

  ! Global vectors
  lat1(:) = 0.0
  do j = 1, ny
     do i = 1, nx
        if (gphit(i,j) > 0.00001) then
           lat1(j) = gphit(i,j)
        endif
     enddo
  enddo

  lon1(:) = 0.0
  do j = 1, ny
     do i = 1, nx
        if (abs(glamt(i,j)) > 0.00001) then
           lon1(i) = glamt(i,j)
        endif
     enddo
  enddo

  ! Local vectors
  lat1_p = lat1(nldj: nlej)
  lon1_p = lon1(nldi: nlei)

  ! Dimension of 2D surface box
  dim_2d = jpiglo*jpjglo
  dim_2d_p = ni_p*nj_p

  if (mype==0 .and. screen>0) then
     write (*,'(/a,5x,a)') 'NEMO-PDAF', '*** NEMO: grid dimensions ***' 
     write(*,'(a,3x,2(6x,a),9x,a)') 'NEMO-PDAF', 'jpiglo','jpjglo','jpk' 
     write(*,'(a,3x,3i12)') 'NEMO-PDAF', jpiglo, jpjglo, jpk
     write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Dimension of global 3D grid box', jpiglo*jpjglo*jpk
     write(*,'(a,5x,a,i12)') 'NEMO-PDAF', 'Number of global surface points', dim_2d
  end if

  if (npes>1 .and. screen>1) then
     write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
          'NEMO-PDAF', 'PE', mype, 'Dimension of local 3D grid box', ni_p*nj_p*nz
     write(*,'(a,2x,a,1x,i4,1x,a,i12)') &
          'NEMO-PDAF', 'PE', mype, 'Number of local surface points', dim_2d_p
  end if


! ***************************************************
! *** Initialize index arrays for wet grid points ***
! ***************************************************

  ! Count number of wet points
  cnt = 0
  do k = 1, nz
     do j = 1, nj_p
        do i = 1, ni_p 
           cnt = cnt + 1
           if (abs(tmp_4d(i,j,k,1)- missing_value) > 0.1) then
              if (k==1) nwet = nwet + 1
              nwet3d = nwet3d + 1
           endif
        enddo
     enddo
  enddo

  ! Initialize index arrays
  ! - for mapping from nx*ny grid to vector of wet points
  ! - mask for wet points

  allocate(wet_pts(7, nwet))
  call memcount(1, 'r', 5*nwet)

  cnt = 0
  cnt_all = 0
     do j = 1, nj_p
        do i = 1, ni_p 
        cnt_all = cnt_all + 1
        if (abs(tmp_4d(i,j,1,1)- missing_value) > 0.1) then
           cnt = cnt+1
           wet_pts(1,cnt) = i + nldi - 1     ! Global longitude index 
           wet_pts(2,cnt) = j + nldj - 1     ! Global latitude index
           wet_pts(6,cnt) = i                ! Longitude index in subdomain
           wet_pts(7,cnt) = j                ! Latitidue index in subdomain
           cnt_layers = 1
           do k = 2, nz
              if (abs(tmp_4d(i,j,k,1)- missing_value) > 0.1) then
                 cnt_layers = cnt_layers + 1
              endif
           enddo
           wet_pts(3,cnt) = cnt_layers
           wet_pts(4,cnt) = cnt_all
        endif
     enddo
  enddo

  ! row 5 stores wet_pts index for vertical column storage in state vector
  wet_pts(5,1) = 1
  do i = 2 , nwet
     wet_pts(5,i) = wet_pts(5,i-1) + wet_pts(3,i-1)
  end do

  ! Initialize index arrays 
  ! - for mapping from vector of wet points to 2d box
  ! - for mapping from vector model grid point coordinats to wet point index
  ! - for retrieving number of wet layers at a grid point coordinate
  ! these arrays also serve as mask arrays (0 indicates land point)

  allocate(idx_wet_2d(ni_p, nj_p))
  allocate(idx_nwet(ni_p, nj_p))
  allocate(nlev_wet_2d(ni_p, nj_p))
  call memcount(1, 'i', 3*ni_p*nj_p)

  idx_wet_2d = 0
  idx_nwet = 0
  nlev_wet_2d = 0
  do i = 1 , nwet
!     idx_wet_2d(wet_pts(1,i), wet_pts(2,i)) = wet_pts(4,i)
!     nlev_wet_2d(wet_pts(1,i), wet_pts(2,i)) = wet_pts(3,i)
     idx_wet_2d(wet_pts(6,i), wet_pts(7,i)) = wet_pts(4,i)
     nlev_wet_2d(wet_pts(6,i), wet_pts(7,i)) = wet_pts(3,i)
  end do
  if (use_wet_state/=2) then
     do i = 1 , nwet
        idx_nwet(wet_pts(6,i), wet_pts(7,i)) = i
     end do
  else
     do i = 1 , nwet
        idx_nwet(wet_pts(6,i), wet_pts(7,i)) = wet_pts(5,i)
     end do
  end if

  if (npes==1) then
     write(*,'(a,5x,a,3x,i11)') 'NEMO-PDAF', 'Number of wet surface points', nwet
     write(*,'(a,5x,a,8x,i11)') 'NEMO-PDAF', 'Number of 3D wet points', nwet3d
     write(*,'(a,5x,a,8x,i11)') 'NEMO-PDAF', '2D wet points * nlayers', nwet*nz
  else 
     if (screen>1) then
        write(*,'(a,2x,a,1x,i4,2x,a,3x,i11)') &
             'NEMO-PDAF', 'PE', mype, 'Number of wet surface points', nwet
        write(*,'(a,2x,a,1x,i4,2x,a,8x,i11)') &
             'NEMO-PDAF', 'PE', mype, 'Number of 3D wet points', nwet3d
        write(*,'(a,2x,a,1x,i4,2x,a,8x,i11)') &
             'NEMO-PDAF', 'PE', mype, '2D wet points * nlayers', nwet*nz
     end if
  end if

  call MPI_Barrier(MPI_COMM_WORLD, MPIerr)


! ********************************************************
! *** Set dimension of 2D and 3D field in state vector ***
! ********************************************************
  
    if (use_wet_state==1) then
       ! State vector contains full columns when surface grid point is wet
       sdim3d = nwet*jpk
       sdim2d = nwet
    elseif (use_wet_state==2) then
       ! State vector only contains wet grid points
       sdim3d = nwet3d
       sdim2d = nwet
    else
       ! State vector contains 3d grid box
       sdim3d = jpiglo*jpjglo*jpk
       sdim2d = jpiglo*jpjglo
    end if

end subroutine initialize
