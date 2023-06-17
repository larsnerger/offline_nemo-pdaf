!$Id: obs_prof_pdafomi.F90 561 2020-11-21 09:59:46Z lnerger $
!> PDAF-OMI observation module for type prof observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!!
!! __Observation type prof:__
!! Handle profile observations (prepared for SHARK data)
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
module obs_prof_pdafomi

  use mod_parallel_pdaf, &
       only: mype=>mype_filter, npes=>npes_filter    ! Rank of filter process
  use PDAFomi, &
       only: obs_f, obs_l   ! Declaration of observation data types
 
  implicit none
  save

  ! Variables which are inputs to the module (usually set in init_pdaf)
  logical :: assim_prof        !< Whether to assimilate this data type
  real    :: rms_obs_prof      !< Observation error standard deviation (for constant errors)
  real(8) :: lradius_prof      !< Localization radius
  real(8) :: sradius_prof      !< Support radius for weight function
  character(len=200)  :: path_prof        ! File path of observation
  character(len=80)   :: file_prof        ! File name of observation
  character (len = 20) :: variable_prof   ! Name of variable to be assimilated

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
  subroutine init_dim_obs_prof(step, dim_obs)

    use PDAFomi, &
         only: PDAFomi_gather_obs
    use mod_assimilation_pdaf, &
         only: filtertype, &
         shiftObsInWet, GaussTransf, trafoConst
    use mod_statevector_pdaf, &
         only: sfields, id
    use mod_nemo_pdaf, only: lat1, lon1, nlats=>jpjglo, nlons=>jpiglo, nlvls=>jpk, &
         lons=>glamt, lats=>gphit, depths=>gdept_1d, idx_wet_2d, &
         nlev_wet_2d, idx_nwet, nwet, dim_2d, use_wet_state, &
         nldi, nldj, nlei, nlej, dim_2d_p, ni_p

    implicit none

! *** Arguments ***
    integer, intent(in)    :: step          !< Current time step
    integer, intent(inout) :: dim_obs       !< Dimension of full observation vector

! *** Local variables ***
    integer(4) :: i, j, k                   ! Counters
    integer(4) :: cnt, cnt0                 ! Counters
    integer(4) :: dim_obs_p                 ! Number of process-local observations
    real(8), allocatable :: obs_p(:)        ! PE-local observation vector
    real(8), allocatable :: ivar_obs_p(:)   ! PE-local inverse observation error variance
    real(8), allocatable :: ocoord_p(:,:)   ! PE-local observation coordinates 
    integer(4) ::  ios                      ! iostat for file opening
    character(len=250) :: obsfile           ! String for observation file name with the path
    real(8)    :: lonp, latp, depthp, obsp  ! Values read from file
    integer(4) :: id_lat, id_lon, id_depth  ! Indices of obs. in model 3D grid box
    real(8)    :: dlat_min, dlon_min, ddepth_min  ! Minimum differences
    integer(4) :: d, i2, j2, k2             ! Counters for shifting
    integer(4) :: shift_i, shift_j, shift_k ! Indices after shifting  
    logical :: wetpoints                    ! Flag for wet points
    real(8) :: minDist, dist                ! distances in shift algorithm
    logical :: debug = .false.              ! Whether to output debug values
    real(8) :: diff                         ! Index difference
    integer(4) :: id_lon_p, id_lat_p        ! PE-local coordinate indices
    integer(4), allocatable :: obs_wet(:)


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    if (mype==0) &
         write (*,'(8x,a)') 'Assimilate observations - obs type PROFILES'

    ! Store whether to assimilate this observation type (used in routines below)
    if (assim_prof) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 0   ! 0=Cartesian

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2



! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! initialize
    if (shiftObsInWet) then
       open (150, file='shiftObsInWet.dat', status='old')
    endif
  
    ! Count observations
    obsfile = trim(path_prof)//trim(file_prof)
    if (mype==0) write (*,'(8x, a, a)')  'Reading observations from: ', trim(obsfile)
    open(unit=120, file=trim(obsfile),iostat=ios)
    if (ios /= 0) write(*,'(8x, a, a)') 'Could not open file ',obsfile 

    dim_obs_p = 0
    cntobs: do
       read(120,*,iostat=ios) lonp, latp, depthp, obsp
       if (ios/=0) exit

       id_lat = 0
       id_lon = 0
       id_depth = 0
     
       dlat_min = 100.0
       latloop: do j = 1, nlats
          diff = abs(lat1(j)-latp)
          if (diff < dlat_min ) then
             dlat_min = diff 
             id_lat = j
!             exit latloop
          endif
       enddo latloop

       dlon_min = 100.0
       lonloop: do i = 1, nlons
          diff = abs(lon1(i)-lonp)
          if (diff < dlon_min ) then
             dlon_min = diff
             id_lon = i
!             exit lonloop
          endif
       enddo lonloop
!       
!        dlat_min=100.0
!        do j = 1, nlats
!           do i = 1, nlons
!              if (sqrt((lats(i,j)-latp)**2.0+(lons(i,j)-lonp)**2.0)< dlat_min) then
!                 dlat_min=sqrt((lats(i,j)-latp)**2.0+(lons(i,j)-lonp)**2.0)
!                 id_lat = j
!                 id_lon = i 
!              endif
!           enddo
!        enddo

       ddepth_min = 100.0
       deploop: do k = 1,nlvls
          diff = abs(depths(k)-depthp)
          if (diff < ddepth_min ) then
             ddepth_min = diff
             id_depth = k 
             exit deploop
          endif
       enddo deploop

!       if (id_lat/=0 .and. id_lon/=0 .and. id_depth/=0) then
       if (id_lat>=nldj .and. id_lat<=nlej &
            .and. id_lon>=nldi .and. id_lon<=nlei .and. id_depth/=0) then
          dim_obs_p = dim_obs_p + 1
       endif

    enddo cntobs

    close(120)
  
    if (npes>1) then
       write (*,'(a, i4, 2x, a, i6)') 'PE', mype, 'Total number of observations ', dim_obs_p
    else
       write (*,'(8x, a, i6)') 'Total number of observations ', dim_obs_p
    end if
  
    dim_obs = dim_obs_p

    if (debug) then
       write(*,*) "       nlons       nlats         nlvls       " 
       write(*,'(i12,i12,i12,i12)') nlons,nlats,nlvls
    end if

    haveobs: if (dim_obs_p > 0) then

       ! Initialize vector of observations and index array
       allocate(obs_p(dim_obs_p))
       allocate(ocoord_p(3, dim_obs_p))
       allocate(obs_wet(dim_obs_p))
       allocate(thisobs%id_obs_p(1, dim_obs_p))

       open(unit=120, file=trim(obsfile),iostat=ios)
       if (ios /= 0) write(*,*) 'Could not open file ',obsfile 

       cnt = 0
       do
          read(120,*,iostat=ios) lonp,latp,depthp,obsp
          if (ios/=0) exit

          if (debug) write(*,'(a,2f10.4)')'Observation: latp/lonp', latp, lonp

          id_lat = 0
          id_lon = 0
          id_depth = 0
      
          !dlat_min = 100.0
          !do j = 1, nlats
        !  if (abs(lat1(j)-latp) < dlat_min ) then
	!    dlat_min = abs(lat1(j)-latp)
	!    id_lat = j
	!  endif
        !enddo

        !dlon_min = 100.0
        !do i = 1, nlons
	!  if (abs(lon1(i)-lonp) < dlon_min ) then
	!    dlon_min = abs(lon1(i)-lonp)
	!    id_lon = i
	!  endif
          !enddo
          dlat_min=10000.0
          do j = 1, nlats
             do i = 1, nlons
                diff = (lats(i,j)-latp)**2.0+(lons(i,j)-lonp)**2.0
                if (diff < dlat_min) then
                   dlat_min=diff
                   id_lat = j
                   id_lon = i 
                endif
             enddo
          enddo

          ddepth_min = 100.0
          do k = 1,nlvls
             diff = abs(depths(k)-depthp)
             if (diff < ddepth_min ) then
                ddepth_min = diff
                id_depth = k 
             endif
          enddo

          if (debug) then
             write(*,'(a,3i6)')'id_lat id_lon id_depth',id_lat,id_lon,id_depth 
             write(*,'(a,3f10.4)')'NEMO: lat lon depth',lat1(id_lat),lon1(id_lon),depths(id_depth)
          end if


       insubdomain: if (id_lat>=nldj .and. id_lat<=nlej &
            .and. id_lon>=nldi .and. id_lon<=nlei .and. id_depth/=0) then

          ! Set i,j,k according to ids - jus to be able to use the existing code
          i = id_lon
          j = id_lat
          k = id_depth

          id_lon_p = id_lon - nldi + 1
          id_lat_p = id_lat - nldj + 1

          ! Compute index of obs. in state vector
          if (use_wet_state==1) then
             cnt0 = idx_nwet(id_lon_p, id_lat_p) + nwet*(k-1)
          elseif (use_wet_state==2) then
             cnt0 = idx_nwet(id_lon_p, id_lat_p) + k - 1
          else
             cnt0 = id_lon_p + ni_p*(id_lat_p-1) + dim_2d_p*(k-1)
          end if

          cnt = cnt + 1

          thisobs%id_obs_p(1, cnt) = cnt0 + sfields(id%temp)%off
          obs_p(cnt) = obsp
          ocoord_p(1,cnt) = real(id_lon)  ! lon
          ocoord_p(2,cnt) = real(id_lat)  ! lat
          ocoord_p(3,cnt) = real(id_depth)

          ! Check whether surface grid point is wet
          wetsurf: if (idx_wet_2d(id_lon_p, id_lat_p)>0) then

             ! If surface is wet check if grid point below surface is wet
             if (id_depth<=nlev_wet_2d(id_lon_p, id_lat_p)) then
                obs_wet(cnt)=1
             else
                obs_wet(cnt)=0
                write(*,*)'+++ Observation point at:',i,j,k,'is not wet point.'
                if (shiftObsInWet) then
                   write(*,*)'Shift Algorithm for observation at point:',i,j,k,'is started.'
                   d=1
                   wetpoints=.false.
                   minDist=100.0D0
                   shift_i=i
                   shift_j=j
                   shift_k=k
                   do while(.not.wetpoints)
                      do k2=max((k-d),1),min((k+d),nlvls)
                         do j2=max((j-d),1),min((j+d),nlats)
                            do i2=max((i-d),1),min((i+d),nlons)
                               wetsurf2: if (idx_wet_2d(i2, j2)>0) then
                                  if (k2<=nlev_wet_2d(i2, j2)) then

!                               if (wet_points_p((k2-1)*nlons*nlats+(j2-1)*nlons+i2)==1)then
                                     write(*,*)'Wet point found at:',i2,j2,k2
                                     wetpoints=.true.
                                     dist=sqrt((lats(i2,j2)-latp)**2.0D0+(lons(i2,j2)-lonp)**2.0D0+(depths(k2)-depthp)**2.0D0)
                                     write(*,*)'Distance to observation point is',dist
                                     if (dist< minDist)then
                                        write(*,*)'Minimum distance is reset to',dist
                                        minDist=dist
                                        shift_i=i2
                                        shift_j=j2
                                        shift_k=k2
                                     endif
                                  endif
                               end if wetsurf2
                            enddo
                         enddo
                      enddo
                      d=d+1
                   enddo
                !here shift
                   write(*,*)'Observation is shifted from',i,j,k,'to',shift_i,shift_j,shift_k 
                !write(150,*)i,j,k,shift_i,shift_j,shift_k,lats(i,j),lons(i,j),depths(k),lats(shift_i,shift_j)&
!           ,lons(shift_i,shift_j),depths(shift_k)
                   write(150,*)i,j,k,shift_i,shift_j,shift_k
                   ocoord_p(1,cnt) = real(shift_i)
                   ocoord_p(2,cnt) = real(shift_j)
                   ocoord_p(3,cnt) = real(shift_k)
                   obs_wet(cnt)=1

                   id_lon_p = shift_i - nldi + 1
                   id_lat_p = shift_j - nldj + 1

                   ! Compute index of obs. in state vector
                   if (use_wet_state==1) then
                     cnt0 = idx_nwet(id_lon_p, id_lat_p) + nwet*(shift_k-1)
                   elseif (use_wet_state==2) then
                     cnt0 = idx_nwet(id_lon_p, id_lat_p) + shift_k - 1
                   else
                     cnt0 = id_lon_p + ni_p*(id_lat_p-1) + dim_2d_p*(shift_k-1)
                   end if

                   thisobs%id_obs_p(1, cnt) = cnt0 + sfields(id%temp)%off
                endif !shiftObsInWet
          
                if (debug) then
                   write(*,'(a,3i6)')'shift: id_lat id_lon id_depth',id_lat,id_lon,id_depth 
                   write(*,'(a,3f10.3)')'shift: lat lon depth nemo',lat1(id_lat),lon1(id_lon),depths(id_depth)
                   write(*,'(a,3f10.3)') 'ocoord_p:',ocoord_p(1,cnt),ocoord_p(2,cnt),ocoord_p(3,cnt)
                end if
             endif ! Wet 3d
          end if wetsurf
       end if insubdomain
       enddo ! read the file

       close(120)

       if (shiftObsInWet) close(150)
    
       if (variable_prof == 'OXY') then
          obs_p = obs_p*44.66
       endif

       select case (GaussTransf)
       case(0)
          write(*,*) 'No Transformation of bio variable'
       case(1)
          write(*,*) 'use log basis 10 transformation of bio variable'
          obs_p = log10(obs_p + trafoConst)
       case(2)
          write(*,*) 'use ln transformation of bio variable'
          obs_p = log(obs_p + trafoConst)
       case(3)
          write(*,*) 'no transformation- box cox still needs to be implemented'
       case DEFAULT
          write(*,*) 'No Transformation of bio variable'
       end select

    else haveobs

       ! If there are no observations just allocate arrays at minimum size
       allocate(obs_p(1))
       allocate(ocoord_p(3, 1))
       allocate(obs_wet(1))
       allocate(thisobs%id_obs_p(1, 1))

    endif haveobs


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ! *** Set inverse observation error variances ***

    allocate(ivar_obs_p(dim_obs_p))

    ivar_obs_p(:) = 1.0 / (rms_obs_prof*rms_obs_prof)


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    call PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         3, lradius_prof, dim_obs)


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

!     IF (twin_experiment .AND. filtertype/=100) THEN
!        CALL read_syn_obs(file_syntobs_TYPE, dim_obs, thisobs%obs_f, 0, 1-mype)
!     END IF


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    deallocate(obs_p, ocoord_p, ivar_obs_p)

  end subroutine init_dim_obs_prof



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
  subroutine obs_op_prof(dim_p, dim_obs, state_p, ostate)

    use PDAFomi, &
         only: PDAFomi_obs_op_gridpoint

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
       ! observation operator for observed grid point values
       call PDAFomi_obs_op_gridpoint(thisobs, state_p, ostate)
    end if

  end subroutine obs_op_prof



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
  subroutine init_dim_obs_l_prof(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    use PDAFomi, only: PDAFomi_init_dim_obs_l

    ! Include localization radius and local coordinates
    use mod_assimilation_pdaf, &
         only: domain_coords, locweight

    implicit none

! *** Arguments ***
    integer, intent(in)  :: domain_p     !< Index of current local analysis domain
    integer, intent(in)  :: step         !< Current time step
    integer, intent(in)  :: dim_obs      !< Full dimension of observation vector
    integer, intent(inout) :: dim_obs_l  !< Local dimension of observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    call PDAFomi_init_dim_obs_l(thisobs_l, thisobs, domain_coords, &
         locweight, lradius_prof, sradius_prof, dim_obs_l)

  end subroutine init_dim_obs_l_prof



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
  subroutine localize_covar_prof(dim_p, dim_obs, HP_p, HPH, coords_p)

    ! Include PDAFomi function
    use PDAFomi, only: PDAFomi_localize_covar

    ! Include localization radius and local coordinates
!    USE mod_assimilation, &   
!         ONLY: locweight

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

!    CALL PDAFomi_localize_covar(thisobs, dim_p, locweight, lradius_prof, srange, &
!         coords_p, HP_p, HPH)

  end subroutine localize_covar_prof

end module obs_prof_pdafomi
