!$Id: obs_TYPE_pdafomi_TEMPLATE.F90 579 2020-11-23 07:32:00Z lnerger $
!> PDAF-OMI template observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
!!
!! Observation type: SST observations from CMEMS
!! 
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF.
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs)
!! and for the observation operator (obs_op).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!!
!! **Using this template:**
!!   To be able to distinguish the observation type and the routines in this module,
!!   we recommend to rename the module according to the observation module-type.
!!   Further,we recommend to replace 'TYPE' in the routine names according to the
!!   type of the observation so that they can be identified when calling them from 
!!   the call-back routines.
!!
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! In addition, there are two optional routine, which are required if filters 
!! with localization are used:
!! * init_dim_obs_l_TYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!! * localize_covar_TYPE \n
!!           Only required if the localized EnKF is used:
!!           Apply covariance localization in the LEnKF.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
module obs_sst_cmems_pdafomi

  use mod_parallel_pdaf, &
       only: mype_filter    ! Rank of filter process
  use PDAFomi, &
       only: obs_f, obs_l   ! Declaration of observation data types
 
  implicit none
  save

  ! Variables which are inputs to the module (usually set in init_pdaf)
  logical :: assim_sst_cmems=.false. !< Whether to assimilate this data sst_cmems
  real    :: rms_obs_sst_cmems       !< Observation error standard deviation (for constant errors)
  real(8) :: lradius_sst_cmems       !< Localization radius
  real(8) :: sradius_sst_cmems       !< Support radius for weight function
  integer :: mode_sst_cmems=0        !< Observation mode: 
                                     !< (0) linear interpolation
                                     !< (1) super-obbing: average 4 observation values
  integer :: time_sst_cmems=0        !< Time at which the daily observations are assimilated
                                     !< (0) mightnight, (12) noon
  character(len=3) :: dist_sst_cmems = 'geo'  ! Type of distance computation: 
                                     !< (gp) for Cartesian distance in unit of grid points
                                     !< (geo) for geographic distance in km
  character(len=200) :: path_sst_cmems  !< Path to CMEMS SST data
  character(len=200) :: file_sst_cmems  !< Filename of CMEMS SST data
  
  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.


! ***********************************************************************
! *** The following two data types are used in PDAFomi                ***
! *** They are declared in PDAFomi and only listed here for reference ***
! ***********************************************************************

! Data type to define the full observations by internally shared variables of the module
!   TYPE obs_f
!           Mandatory variables to be set in INIT_DIM_OBS
!      INTEGER :: doassim                   ! Whether to assimilate this observation type
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!                                           ! (0) Cartesian, (1) Cartesian periodic
!                                           ! (2) simplified geographic, (3) geographic haversine function
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! Indices of observed field in state vector (process-local)
!           
!           Optional variables - they can be set in INIT_DIM_OBS
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!
!           Variables with predefined values - they can be changed in INIT_DIM_OBS
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!      INTEGER :: use_global_obs=1          ! Whether to use (1) global full obs. 
!                                           ! or (0) obs. restricted to those relevant for a process domain
!
!           The following variables are set in the routine PDAFomi_gather_obs
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      INTEGER :: dim_obs_g                 ! global number of observations
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
!      INTEGER :: off_obs_g                 ! Offset of this observation in overall global obs. vector
!      INTEGER :: obsid                     ! Index of observation over all assimilated observations
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER, ALLOCATABLE :: id_obs_f_lim(:) ! Indices of domain-relevant full obs. in global vector of obs.
!                                           ! (only if full obs. are restricted to process domain))
!   END TYPE obs_f

! Data type to define the local observations by internally shared variables of the module
!   TYPE obs_l
!      INTEGER :: dim_obs_l                 ! number of local observations
!      INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
!      REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
!      REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
!      INTEGER :: locweight                 ! Specify localization function
!      REAL :: lradius                      ! localization radius
!      REAL :: sradius                      ! support radius for localization function
!   END TYPE obs_l
! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could renamed the variables
  type(obs_f), target, public :: thisobs      ! full observation
  type(obs_l), target, public :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

contains

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
!! **Adapting the template**
!! In this routine the variables listed above have to be initialized. One
!! can include modules from the model with 'use', e.g. for mesh information.
!! Alternatively one could include these as subroutine arguments
!!
  subroutine init_dim_obs_sst_cmems(step, dim_obs)

    use PDAFomi, &
         only: PDAFomi_gather_obs, PDAFomi_get_interp_coeff_lin
    use mod_assimilation_pdaf, &
         only: filtertype, iday, deg2rad
    use mod_statevector_pdaf, &
         only: id, sfields
    use mod_parallel_pdaf, only: mype_filter, npes_filter
    use mod_io_pdaf, only: check
    use mod_nemo_pdaf, only: lat1_p, lon1_p, nlats=>nj_p, nlons=>ni_p, &
         idx_nwet, use_wet_state, nlei, nlej
    use netcdf

    implicit none

! *** Arguments ***
    integer, intent(in)    :: step       !< Current time step
    integer, intent(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    logical :: debug = .false.           ! Activate debugging output for index calculations
    integer :: i, j, cnt                 ! Counters
    integer :: ido_start, ido_end        ! Counters
    integer :: idm_start, idm_end        ! Counters
    integer :: dim_obs_p                 ! Number of process-local observations
    real, allocatable :: obs_p(:)        ! PE-local observation vector
    real, allocatable :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    real, allocatable :: ocoord_p(:,:)   ! PE-local observation coordinates 
    logical :: doassim_now=.true.        ! Whether we assimilate the observation at the current time
    integer(4) :: status                 ! Status flag for availability of observations
    character(len=100) :: file_full      ! filename including path
    character(len=2) :: strday           ! day as string
    integer(4) :: ncid, dimid, lonid, latid, varid  ! nc file IDs
    integer(4) :: startv(3), cntv(3)                ! Index arrays for reading from nc file
    integer(4) :: dim_olat, dim_olon                ! Grid dimensions read from file
    integer(4), allocatable :: sst_aux(:,:)         ! SST field read from file 
    real, allocatable :: lon_obs(:), lat_obs(:)     ! Obs. coordinates read from file
    real, allocatable :: lon_model(:), lat_model(:) ! Longitude/latitude of model in radians
    integer :: iderr                                ! Error flag for determining indices
    real, parameter :: sst_scale = 0.01             ! Scaling factor to convert file value to degC
    integer(4) :: ido_n, ido_e, ido_s, ido_w        ! Obs. ID limits N/E/S/W for model grid
    integer :: idm_n, idm_e, idm_s, idm_w           ! Model ID limits NESW 
    real :: wlonM, elonM, nlatM, slatM              ! Coordinate limits of model grid
    real :: dlonM, dlatM                            ! Model grid spacing
    real :: dlonO, dlatO                            ! Observation grid spacing
    real :: latM_limit                              ! Comparison limit in latitude for model coordinate
    real :: lonM, latM                              ! Longitude/latitude of a model grid point
    real :: gcoords(4,2)                  ! Grid point coordinates for computing interpolation coeffs
    integer(4) :: obsflag                 ! Count observation in direct vicinity
    integer(4) :: cntobs(4)               ! Count grid points with 0 to 4 obs. neighbours
    integer(4) :: obs_sum                 ! Sum of observation integer values
    integer :: sgn_olat                   ! Orientation of latitude Obs.: -1 for north-south/+1 for south-north
    integer :: sgn_mlat                   ! Orientation of latitude Model: -1 for north-south/+1 for south-north


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    if (mype_filter==0) &
         write (6,'(8x,a)') 'Assimilate observations - OBS_SST_CMEMS'

    ! Store whether to assimilate this observation type (used in routines below)
    if (assim_sst_cmems) thisobs%doassim = 1

    ! Specify type of distance computation
    if (trim(dist_sst_cmems) == 'gp') then
       if (mype_filter==0) write (*,'(8x,a)') '--- use Cartesian grid point distances'
       thisobs%disttype = 0   ! 0=Cartesian
    else
       if (mype_filter==0) write (*,'(8x,a)') '--- use geographic distances'
       thisobs%disttype = 2   ! 2=Geographic
    end if

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! In case of MPI parallelization restrict observations to sub-domains
    if (npes_filter>1) thisobs%use_global_obs = 0


! **********************************
! *** Read PE-local observations ***
! **********************************

    doassim: if (doassim_now) then

       ! Only execute this if we assimilate observations at this hour

       ! The SST data can can be downloaded as follows:
       ! python -m motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu 
       !    --service-id SST_EUR_SST_L3S_NRT_OBSERVATIONS_010_009_a-TDS 
       !    --product-id METEOFRANCE-EUR-SST_L3MULTISENSOR_NRT-OBS_FULL_TIME_SERIE 
       !    --longitude-min -4.5 --longitude-max 30.5 --latitude-min 48.5 --latitude-max 66 
       !    --date-min "2018-10-01 00:00:00" --date-max "2018-10-31 00:00:00" 
       !    --variable adjusted_sea_surface_temperature 
       !    --out-name sst_multi_201810.nc --user <USERNAME> --pwd <PASSWD>


       ! read observation values and their coordinates
       file_full = trim(path_sst_cmems)//trim(file_sst_cmems)

       if (mype_filter==0) write (6,'(8x,a,i3,a)') 'Read observations for day ',iday,' from file:'
       if (mype_filter==0) write (6,'(10x,a)') trim(file_full)

       ! Open the file. NF90_NOWRITE tells netCDF to have read-only access to file.
       call check( nf90_open(file_full, NF90_NOWRITE, ncid) )

       ! Read dimensions of observation grid
       call check( nf90_inq_dimid(ncid, "lat", dimid) )
       call check( nf90_Inquire_dimension(ncid, dimid, len=dim_olat) )
       call check( nf90_inq_dimid(ncid, "lon", dimid) )
       call check( nf90_Inquire_dimension(ncid, dimid, len=dim_olon) )

       ! Allocate arrays
       allocate(sst_aux(dim_olon, dim_olat))
       allocate(lon_obs(dim_olon), lat_obs(dim_olat))

       ! Get variable IDs and read data
       call check( nf90_inq_varid(ncid, "adjusted_sea_surface_temperature", varid) )
       call check( nf90_inq_varid(ncid, "lon", lonid) )
       call check( nf90_inq_varid(ncid, "lat", latid) )

       call check( nf90_get_var(ncid, lonid, lon_obs) )
       call check( nf90_get_var(ncid, latid, lat_obs) )

       ! Read SST values
       ! They are in deg C but have to be scaled by sst_scale (1/100).
       startv(1) = 1 ! lon
       startv(2) = 1 ! lat
       startv(3) = iday ! time
       cntv(1) = dim_olon
       cntv(2) = dim_olat
       cntv(3) = 1
       call check( nf90_get_var(ncid, varid, sst_aux, start=startv, count=cntv) )

       ! Close the file
       call check( nf90_close(ncid) )

  
       ! *** Convert observation coordinates to radians ***
       lon_obs = lon_obs * deg2rad
       lat_obs = lat_obs * deg2rad

       ! *** Store model coordinates in radians ***
       allocate(lon_model(nlons))
       allocate(lat_model(nlats))
       lon_model = lon1_p * deg2rad
       lat_model = lat1_p * deg2rad

       cartdist: if (trim(dist_sst_cmems) == 'gp') then

          ! Set grid point coordinates for Cartesian distance computation

          dlatM = lat_model(2) - lat_model(1)  ! Model grid spacing in latitude
          dlonM = lon_model(2) - lon_model(1)  ! Model grid spacing in longitude

          do i = 1, dim_olon
             lon_obs(i) = (lon_obs(i) - lon_model(1))/dlonM + 1.0
          end do
          do j = 1, dim_olat
             lat_obs(j) = (lat_obs(j) - lat_model(1))/dlatM + 1.0
          end do

          do i = 1, nlons
             lon_model(i) = i + nlei - 1
          end do
          do j = 1, nlats
             lat_model(j) = j + nlej - 1
          end do

       end if cartdist


! *******************************************************************
! *** Determine indices of observation grid inside the model grid ***
! *******************************************************************

       ! Set index limits for model grid
       idm_n = nlats
       idm_s = 1
       idm_w = 1
       idm_e = nlons

       ! boundaries and spacing for model grid
       nlatM = lat_model(idm_n)             ! Northern latitude limit of the grid box
       slatM = lat_model(idm_s)             ! Southern latitude limit of the grid box
       wlonM = lon_model(idm_w)             ! Western longitude limit of the grid box 
       elonM = lon_model(idm_e)             ! Eastern longitude limit of the grid box
       dlatM = lat_model(2) - lat_model(1)  ! Model grid spacing in latitude
       dlonM = lon_model(2) - lon_model(1)  ! Model grid spacing in longitude
       sgn_mlat = int(sign(1.0,dlatM))      ! Orientation of latitudinal direction (-1: north-to-south)

       ! Observation grid spacing and oriantation
       dlatO = lat_obs(2) - lat_obs(1)      ! Observation grid spacing in latitude
       dlonO = lon_obs(2) - lon_obs(1)      ! Observation grid spacing in longitude
       sgn_olat = int(sign(1.0,dlatO))      ! Orientation of latitudinal direction (-1: north-to-south)

       if (debug) then
          write (*,*) 'dim_olat/lon', dim_olat, dim_olon
          write (*, *) 'idm all', idm_w, idm_e, idm_n, idm_s
          write (*,'(a,4f12.5)') 'mcoords limits',wlonM, elonM, nlatM, slatM
          if (sgn_mlat<0) then
             write (*,'(a,4f12.5)') 'ocoords limits', lon_obs(1), lon_obs(dim_olon), lat_obs(1), lat_obs(dim_olat)
          else
             write (*,'(a,4f12.5)') 'ocoords limits', lon_obs(1), lon_obs(dim_olon), lat_obs(dim_olat), lat_obs(1)
          end if
          write (*,*) 'sgn_olat/sgn_mlat', sgn_olat, sgn_mlat
          write (*,'(a,2es10.2,a,2es10.2)') 'dlat/lon: M:', dlatM, dlonM, ', O:', dlatO, dlonO
       end if


       ! Initialize error flag
       iderr = 0

       ! Compute obs. coordinate indices for model grid limits 
       ! If needed adapt index limits for model grid
       ido_w = ceiling(abs(wlonM - lon_obs(1)) / dlonO)+1

       if (ido_w<2) then
          write (*,*) 'reset ido_w'
          ido_w = 2
          idm_w = ceiling((lon_obs(ido_w) - wlonM) / dlonM) + 1
       elseif (ido_w>dim_olon) then
          iderr = 1
       endif
       if (mode_sst_cmems==0) then
          idm_w = ceiling((lon_obs(ido_w) - wlonM) / dlonM) + 1
       end if

       ido_e = ceiling((elonM - lon_obs(1)) / dlonO)
       if (ido_e>=dim_olon) then
          write (*,*) 'reset ido_e'
          ido_e = dim_olon-1
          idm_e = floor((lon_obs(dim_olon) - wlonM) / dlonM) + 1
       elseif (ido_e<1) then
          iderr = 2
       endif
       if (mode_sst_cmems==0) then
          idm_e = floor((lon_obs(ido_e) - wlonM) / dlonM) + 1
       end if

       if (sgn_olat > 0) then
          ido_n = ceiling((nlatM - lat_obs(1)) / dlatO)
          if (ido_n>dim_olat) then
             write (*,*) 'reset ido_n'
             ido_n = dim_olat-1
             idm_n = ceiling((lat_obs(ido_n) - nlatM) / dlatM) + 1
          elseif (ido_n<1) then
             iderr = 3
          end if
          if (mode_sst_cmems==0) then
             if (sgn_mlat > 0) then
                idm_n = floor(abs(lat_obs(ido_n) - slatM) / dlatM) + 1
             else
                idm_n = ceiling((lat_obs(ido_n) - nlatM) / dlatM) + 1
             end if
          end if

          ido_s = ceiling(abs(slatM - lat_obs(1)) / dlatO) + 1
          if (ido_s<2) then
             write (*,*) 'reset ido_s'
             ido_s = 2
             if (sgn_mlat > 0) then
                idm_s = floor(abs(lat_obs(ido_s) - slatM) / dlatM) + 1
             else
                idm_s = floor(abs(lat_obs(ido_s) - nlatM) / dlatM) + 1
             end if
          elseif (ido_s>dim_olat) then
             iderr = 4
          end if
          if (mode_sst_cmems==0) then
             if (sgn_mlat > 0) then
                idm_s = ceiling(abs(lat_obs(ido_s) - slatM) / dlatM) + 1
             else
                idm_s = floor((lat_obs(ido_s) - nlatM) / dlatM) + 1
             end if
          end if
       else
          ido_n = ceiling((nlatM - lat_obs(1)) / dlatO)+1
          if (ido_n<2) then 
             write (*,*) 'reset ido_n'
             ido_n = 2
             idm_n = floor((lat_obs(ido_n) - nlatM) / dlatM) + 1
          elseif (ido_n>dim_olat) then
             iderr = 5
          endif
          if (mode_sst_cmems==0) then
             if (sgn_mlat > 0) then
                idm_n = floor((lat_obs(ido_n) - slatM) / dlatM) + 1
             else
                idm_n = ceiling((lat_obs(ido_n) - nlatM) / dlatM) + 1
             end if
          end if

          ido_s = ceiling((slatM - lat_obs(1)) / dlatO)
          if (ido_s>dim_olat) then
             write (*,*) 'reset ido_s'
             ido_s = dim_olat-1
             idm_s = ceiling((lat_obs(dim_olat) - nlatM) / dlatM)+1
          elseif (ido_s<1) then
             iderr = 6
          endif
          if (mode_sst_cmems==0) then
             if (sgn_mlat > 0) then
                idm_s = floor(abs(lat_obs(ido_s) - slatM) / dlatM)
             else
                idm_s = floor((lat_obs(ido_s) - nlatM) / dlatM)+1
             end if
          end if
       end if
       if (iderr/=0) write (*,*) 'ERROR: ',iderr,'Observations not overlapping with model grid'

       if (debug) then
          write (*,'(a,2x,4i10)')   'ido in WENS', ido_w, ido_e, ido_n, ido_s
          write (*,'(a,4f12.5)') 'ocoords WENS in ', &
               lon_obs(ido_w), lon_obs(ido_e), lat_obs(ido_n), lat_obs(ido_s)
          write (*,'(a,4f12.5)') 'ocoords WENS out', &
               lon_obs(ido_w-1), lon_obs(ido_e+1), lat_obs(ido_n+sgn_olat), lat_obs(ido_s-sgn_olat)
          
          write (*,'(a,2x,4i10)') 'idm new limits', idm_w, idm_e, idm_n, idm_s
          if (mode_sst_cmems==0) &
               write (*,'(a,4f12.5)') 'mcoords out     ', lon_model(idm_w-1), &
               lon_model(idm_e+1), lat_model(idm_n+sgn_mlat), lat_model(idm_s-sgn_mlat)
          write (*,'(a,4f12.5)') 'mcoords in      ', lon_model(idm_w), &
               lon_model(idm_e), lat_model(idm_n), lat_model(idm_s)
          if (sgn_mlat>0) then
             write (*,*) 'limits SN', 1, nlats, 'new', idm_s, idm_n
          else
             write (*,*) 'limits NS', 1, nlats, 'new', idm_n, idm_s
          end if
          write (*,*) 'limits WE', 1, nlons, 'new', idm_w, idm_e
       end if


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

       obsmodeA: if (mode_sst_cmems==0) then

          ! *** Linear interpolation ***

          if (mype_filter==0) &
               write (6,'(8x,a)') '--- use observations with linear interpolation'

          ! *** Count valid observations that lie within the model grid ***
          cnt = 0

          ! Set start and end index for latitude
          if (sgn_olat == 1) then
             ido_start = ido_s
             ido_end = ido_n
          else
             ido_start = ido_n
             ido_end = ido_s
          end if

          ! Set limit index for computing model latitude index
          if (sgn_mlat>0) then
             latM_limit = slatM
          else
             latM_limit = nlatM
          end if

          ! Loop through observation grid
          do j = ido_start, ido_end
             do i = ido_w, ido_e
                if (sst_aux(i,j) > -10000) then

                   ! find model grid point indices corresponding to observation point
                   idm_w = floor((lon_obs(i) - wlonM) / dlonM) + 1
                   idm_e = idm_w + 1
                   idm_n = ceiling(abs(lat_obs(j) - latM_limit) / dlatM) + 1
                   idm_s = idm_n - sgn_mlat

                   ! Check whether all of the neighboring grid points are wet - then use the observation
                   if (idx_nwet(idm_w, idm_n)>0 .and. idx_nwet(idm_e, idm_n)>0 .and. &
                        idx_nwet(idm_w, idm_s)>0 .and. idx_nwet(idm_e, idm_s) > 0) then
                      cnt = cnt + 1
                   end if

                endif
             end do
          end do

       else obsmodeA

          ! *** Super-Obbing ***

          if (mype_filter==0) &
               write (6,'(8x,a)') '--- use observations with super-obbing'

          ! *** Count valid observations that lie within the grid ***
          cnt = 0
          cntobs = 0

          ! Set start and end index for latitude
          if (sgn_olat == 1) then
             idm_start = idm_s
             idm_end = idm_n
          else
             idm_start = idm_n
             idm_end = idm_s
          end if

          ! Loop through model grid
          do j = idm_start, idm_end
             do i = idm_w, idm_e

                ! Model grid point coordinates
                latM = lat_model(j)
                lonM = lon_model(i)
             
                ! Compute observation grid point indices
                ido_w = ceiling(abs(lonM - lon_obs(1)) / dlonO)
                ido_e = ceiling(abs(lonM - lon_obs(1)) / dlonO)+1
                ido_n = ceiling(abs(latM - lat_obs(1)) / dlatO)+1
                ido_s = ceiling(abs(latM - lat_obs(1)) / dlatO)
             
                ! For wet grid points: Check whether there are valid
                ! observations at the observation grid points around the point
                if (idx_nwet(i, j)>0.0) then

                   if (ido_w<1 .or. ido_w>dim_olon .or. ido_e<1 .or. ido_e>dim_olon &
                        .or. ido_n<1 .or. ido_n>dim_olat .or. ido_s<1 .or. ido_s>dim_olat) &
                        write (*,*) 'ido out of range', ido_w, ido_e, ido_n, ido_s

                   obsflag = 0
                   if (sst_aux(ido_w,ido_s)>-10000) obsflag = obsflag + 1
                   if (sst_aux(ido_e,ido_s)>-10000) obsflag = obsflag + 1
                   if (sst_aux(ido_w,ido_n)>-10000) obsflag = obsflag + 1
                   if (sst_aux(ido_e,ido_n)>-10000) obsflag = obsflag + 1

                   ! Count how many of the neighboring observations are valid
                   if (obsflag==4) cntobs(4) = cntobs(4) + 1
                   if (obsflag==3) cntobs(3) = cntobs(3) + 1
                   if (obsflag==2) cntobs(2) = cntobs(2) + 1
                   if (obsflag==1) cntobs(1) = cntobs(1) + 1

                   ! Count grid points with observations in direct vicinity
                   if (obsflag>0) cnt = cnt+1
                end if

             end do
          end do
          write (6, '(12x, a,4i8)') 'grid points with 1/2/3/4 neighbour obs.', cntobs(1:4)

       end if obsmodeA

       ! Set observation dimension
       dim_obs_p = cnt
       dim_obs = cnt 

       if (npes_filter==1) then
          write (6,'(8x, a, i7)') '--- number of observations from SST_CMEMS: ', dim_obs
       else
          write (6,'(a, i4, 2x, a, i7)') 'PE', mype_filter, &
               '--- number of observations from SST_CMEMS: ', dim_obs
       end if

    else doassim

       ! This is for the case that we do not assimilate this data at the current hour
       dim_obs = 0
       dim_obs_p = 0

    end if doassim


    ! *** Initialize vector of observations on the process sub-domain ***
    ! *** Initialize coordinate array of observations                 ***
    ! *** Initialize process local index array                        ***
    if (dim_obs_p > 0) then

       allocate(obs_p(dim_obs_p))
       allocate(ivar_obs_p(dim_obs_p))
       allocate(ocoord_p(thisobs%ncoord, dim_obs_p))

       obsmodeB: if (mode_sst_cmems==0) then

          ! *** Linear interpolation ***

          ! Allocate process-local index array
          allocate(thisobs%id_obs_p(4, dim_obs_p))

          ! Allocate array of interpolation coefficients. As ID_OBS_P, the number
          ! of rows corresponds to the number of grid points using the the interpolation
          allocate(thisobs%icoeff_p(4, dim_obs_p))

          cnt = 0

          do j = ido_start, ido_end
             do i = ido_w, ido_e

                ! find grid point indices corresponding to observation point
                idm_w = floor((lon_obs(i) - wlonM) / dlonM) + 1
                idm_e = idm_w + 1
                idm_n = ceiling(abs(lat_obs(j) - latM_limit) / dlatM) + 1
                idm_s = idm_n - sgn_mlat

                haveobs: if (sst_aux(i,j) > -10000) then

                   ! find grid point indices corresponding to observation point
                   idm_w = floor((lon_obs(i) - wlonM) / dlonM) + 1
                   idm_e = idm_w + 1
                   idm_n = ceiling(abs(lat_obs(j) - latM_limit) / dlatM) + 1
                   idm_s = idm_n - sgn_mlat

                   ! Use obs. if all of the neighboring grid points are wet
                   wetpoint: if (idx_nwet(idm_w, idm_n)>0 .and. idx_nwet(idm_e, idm_n)>0 .and. &
                        idx_nwet(idm_w, idm_s)>0 .and. idx_nwet(idm_e, idm_s) > 0) then

                      cnt = cnt + 1

                      ! Set indices of 4 grid points neightbors of the observation
                      if (use_wet_state==1 .or. use_wet_state==2) then
                         thisobs%id_obs_p(1, cnt) = idx_nwet(idm_w,idm_s) + sfields(id%temp)%off
                         thisobs%id_obs_p(2, cnt) = idx_nwet(idm_e,idm_s) + sfields(id%temp)%off
                         thisobs%id_obs_p(3, cnt) = idx_nwet(idm_w,idm_n) + sfields(id%temp)%off
                         thisobs%id_obs_p(4, cnt) = idx_nwet(idm_e,idm_n) + sfields(id%temp)%off
                      else
                         thisobs%id_obs_p(1, cnt) = idm_w + nlons*(idm_s-1) + sfields(id%temp)%off
                         thisobs%id_obs_p(2, cnt) = idm_e + nlons*(idm_s-1) + sfields(id%temp)%off
                         thisobs%id_obs_p(3, cnt) = idm_w + nlons*(idm_n-1) + sfields(id%temp)%off
                         thisobs%id_obs_p(4, cnt) = idm_e + nlons*(idm_n-1) + sfields(id%temp)%off
                      end if

                      ! Store observation value and coordinates
                      obs_p(cnt) = real(sst_aux(i, j)) * sst_scale
                      ocoord_p(1, cnt) = lon_obs(i)
                      ocoord_p(2, cnt) = lat_obs(j)

                      ! *** Determine interpolation coefficients ***

                      ! Determine coordinates of grid points around observation
                      ! Order of coefficients:  (3) ---- (4)          
                      !                          |        |
                      !                         (1) ---- (2)
                      ! Only 4 coordinate values are required for bi-linear interpolation
                      gcoords(1,1) = lon_model(idm_w)
                      gcoords(1,2) = lat_model(idm_s)
                      gcoords(2,1) = lon_model(idm_e)
                      gcoords(3,2) = lat_model(idm_n)

                      ! Compute interpolation coefficients
                      call PDAFomi_get_interp_coeff_lin(4, 2, gcoords, ocoord_p(:, cnt), &
                           thisobs%icoeff_p(:,cnt))

                   end if wetpoint
                endif haveobs
             enddo
          enddo

       else obsmodeB

          ! *** Super-Obbing ***

          ! Allocate process-local index array
          allocate(thisobs%id_obs_p(1, dim_obs_p))

          ! Loop through model grid
          cnt = 0
          do j = idm_start, idm_end
             do i = idm_w, idm_e

                ! Model grid point coordinates
                latM = lat_model(j)
                lonM = lon_model(i)
             
                ! Compute observation grid point indices
                ido_w = ceiling(abs(lonM - lon_obs(1)) / dlonO)
                ido_e = ceiling(abs(lonM - lon_obs(1)) / dlonO)+1
                ido_n = ceiling(abs(latM - lat_obs(1)) / dlatO)+1
                ido_s = ceiling(abs(latM - lat_obs(1)) / dlatO)
             
                ! For wet grid points: Check whether there are valid
                ! observations at the observation grid points around the point
                wetpointB: if (idx_nwet(i, j)>0) then

                   obsflag = 0
                   if (sst_aux(ido_w,ido_s)>-10000) obsflag = obsflag + 1
                   if (sst_aux(ido_e,ido_s)>-10000) obsflag = obsflag + 1
                   if (sst_aux(ido_w,ido_n)>-10000) obsflag = obsflag + 1
                   if (sst_aux(ido_e,ido_n)>-10000) obsflag = obsflag + 1

                   ! Use obs. if at least one observation exist in direct vicinity
                   obsflg: if (obsflag>0) then

                      cnt = cnt+1

                      ! Set index of grid point 
                      if (use_wet_state==1 .or. use_wet_state==2) then
                         thisobs%id_obs_p(1, cnt) = idx_nwet(i, j) + sfields(id%temp)%off
                      else
                         thisobs%id_obs_p(1, cnt) = i + nlons*(j-1) + sfields(id%temp)%off
                      end if

                      ! Store observation coordinates
                      ocoord_p(1, cnt) = lon_model(i)
                      ocoord_p(2, cnt) = lat_model(j)

                      ! Compute observation value by averaging
                      obs_sum = 0
                      if (sst_aux(ido_w,ido_s)>-10000) &
                           obs_sum = obs_sum + sst_aux(ido_w,ido_s)
                      if (sst_aux(ido_e,ido_s)>-10000) &
                           obs_sum = obs_sum + sst_aux(ido_e,ido_s)
                      if (sst_aux(ido_w,ido_n)>-10000) &
                           obs_sum = obs_sum + sst_aux(ido_w,ido_n)
                      if (sst_aux(ido_e,ido_n)>-10000) &
                           obs_sum = obs_sum + sst_aux(ido_e,ido_n)
                      obs_p(cnt) = real(obs_sum) * sst_scale / real(obsflag)

                   end if obsflg

                end if wetpointB

             end do
          end do


       end if obsmodeB

    else

       ! for DIM_OBS_P=0

       allocate(obs_p(1))
       allocate(ivar_obs_p(1))
       allocate(ocoord_p(thisobs%ncoord, 1))
       allocate(thisobs%id_obs_p(1, dim_obs_p))

    end if


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ivar_obs_p(:) = 1.0 / (rms_obs_sst_cmems*rms_obs_sst_cmems)


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    ! This routine is generic for the case that only the observations, 
    ! inverse variances and observation coordinates are gathered

    call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, lradius_sst_cmems, dim_obs)


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

!   IF (twin_experiment .AND. filtertype/=11) THEN
!      CALL read_syn_obs(file_syntobs_TYPE, dim_obs, thisobs%obs_f, 0, 1-mype_filter)
!   END IF


! ********************
! *** Finishing up ***
! ********************

    if (doassim_now) then
       ! Deallocate all local arrays
       deallocate(obs_p, ocoord_p, ivar_obs_p)
       deallocate(sst_aux, lon_obs, lat_obs)
    end if

  end subroutine init_dim_obs_sst_cmems



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  subroutine obs_op_sst_cmems(dim_p, dim_obs, state_p, ostate)

    use PDAFomi, &
         only: PDAFomi_obs_op_interp_lin, PDAFomi_obs_op_gridpoint

    implicit none

! *** Arguments ***
    integer, intent(in) :: dim_p                 !< PE-local state dimension
    integer, intent(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    real, intent(in)    :: state_p(dim_p)        !< PE-local model state
    real, intent(inout) :: ostate(dim_obs)       !< Full observed state



! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    if (thisobs%doassim==1) then
    
       if (mode_sst_cmems==0) then

          ! Observation operator for averaging over grid points
          call PDAFomi_obs_op_interp_lin(thisobs, 4, state_p, ostate)

       else

          ! Observation operator for averaging over grid points
          call PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)

       end if
    end if

  end subroutine obs_op_sst_cmems



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  subroutine init_dim_obs_l_sst_cmems(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    use PDAFomi, only: PDAFomi_init_dim_obs_l
    use mod_nemo_pdaf, only: lat1, lon1

    ! Include localization radius and local coordinates
    ! one can also set observation-specific values for the localization.
    use mod_assimilation_pdaf, &   
         only: domain_coords, locweight, deg2rad

    implicit none

! *** Arguments ***
    integer, intent(in)  :: domain_p     !< Index of current local analysis domain
    integer, intent(in)  :: step         !< Current time step
    integer, intent(in)  :: dim_obs      !< Full dimension of observation vector
    integer, intent(inout) :: dim_obs_l  !< Local dimension of observation vector

    real :: coords(2)

! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Here one has to specify the coordinates of the local analysis domain
    ! (coords_l) and the localization variables, which can be different for
    ! each observation type and can be made dependent on the index DOMAIN_P.
    ! coords_l should be set in the call-back routine init_dim_l.
    ! coords_l is domain_coords in HBM

    if (trim(dist_sst_cmems) == 'gp') then
       coords(1) = domain_coords(1)
       coords(2) = domain_coords(2)
    else
       coords(1) = lon1(int(domain_coords(1)))*deg2rad
       coords(2) = lat1(int(domain_coords(2)))*deg2rad
    end if

    call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords, &
         locweight, lradius_sst_cmems, sradius_sst_cmems, dim_obs_l)

  end subroutine init_dim_obs_l_sst_cmems



!-------------------------------------------------------------------------------
!> Perform covariance localization for local EnKF on the module-type observation
!!
!! The routine is called in the analysis step of the localized
!! EnKF. It has to apply localization to the two matrices
!! HP and HPH of the analysis step for the module-type
!! observation.
!!
!! This routine calls the routine PDAFomi_localize_covar
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type.
!!
  subroutine localize_covar_sst_cmems(dim_p, dim_obs, HP_p, HPH, coords_p)

    ! Include PDAFomi function
    use PDAFomi, only: PDAFomi_localize_covar

    ! Include localization radius and local coordinates
    use mod_assimilation_pdaf, &   
         only: locweight

    implicit none

! *** Arguments ***
    integer, intent(in) :: dim_p                 !< PE-local state dimension
    integer, intent(in) :: dim_obs               !< Dimension of observation vector
    real, intent(inout) :: HP_p(dim_obs, dim_p)  !< PE local part of matrix HP
    real, intent(inout) :: HPH(dim_obs, dim_obs) !< Matrix HPH
    real, intent(in)    :: coords_p(:,:)         !< Coordinates of state vector elements



! *************************************
! *** Apply covariance localization ***
! *************************************

    ! Here one has to specify the three localization variables
    ! which can be different for each observation type.

    call PDAFomi_localize_covar(thisobs, dim_p, locweight, lradius_sst_cmems, sradius_sst_cmems, &
         coords_p, HP_p, HPH)

  end subroutine localize_covar_sst_cmems

end module obs_sst_cmems_pdafomi
