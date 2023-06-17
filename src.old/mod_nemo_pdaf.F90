module mod_nemo_pdaf

  ! *** NEMO model variables
  integer :: jpiglo, jpjglo, jpk        ! Global NEMO grid dimensions
  integer :: nldi, nldj                 ! first inner index in i/j direction of sub-domain
  integer :: nlei, nlej                 ! last inner index in i/j direction of sub-domain
  integer :: nimpp, njmpp               ! start i,j of subdomain including halo
  integer :: jptra = 16                 ! Number of prognostic ERGOM tracers
  integer :: jptra2 = 3                 ! Number of diagnostic ERGOM tracers

  real, allocatable :: glamt(:,:)       ! Longitudes
  real, allocatable :: gphit(:,:)       ! Latitudes
  real, allocatable :: gdept_1d(:)      ! Depths (in files `nav_lev`)
  
  ! *** Other grid variables
  real, allocatable :: tmp_4d(:,:,:,:)     ! 4D array used to represent full NEMO grid box
  real, allocatable :: lat1(:), lon1(:)    ! Vectors holding latitude and latitude

  integer :: dim_2d                        ! Dimension of 2d grid box    
  integer :: nwet                          ! Number of surface wet grid points
  integer :: nwet3d                        ! Number of 3d wet grid points
  integer :: sdim2d, sdim3d                ! Dimension of 2D/3D field in state vector
  integer, allocatable :: wet_pts(:,:)     ! Index array for wet grid points
                           ! (1) latitude, (2) langitude, (3) number wet layers, (4) index in 2d grid box
  integer, allocatable :: idx_wet_2d(:,:)  ! Index array for wet_pts row index in 2d box
  integer, allocatable :: idx_nwet(:,:)    ! Index array for wet_pts row index in wet surface grid points
  integer, allocatable :: nlev_wet_2d(:,:) ! Number of wet layers for ij position in 2d box

  integer :: use_wet_state=1

  integer :: i0, j0                         ! PE-local halo offsets
  integer :: istart, jstart                 ! Start indices for internal local domain
  integer :: dim_2d_p, dim_3d_p             ! Dimension of 2d/3d grid box of sub-domain   
  integer :: ni_p, nj_p, nk_p               ! Size of decomposed grid
  real, allocatable :: lat1_p(:), lon1_p(:) ! Vectors holding latitude and latitude for decomposition

  ! *** File name and path to read grid information
  character(len=200)  :: path_dims         ! Path for NEMO file holding dimensions
  character(len=80)   :: file_dims         ! File name NEMO file holding dimensions


end module mod_nemo_pdaf
