!$Id$
!BOP
!
! !ROUTINE: assimilation_pdaf_offline - Control PDAF offline analysis
!
! !INTERFACE:
SUBROUTINE assimilation_pdaf_offline()

! !DESCRIPTION:
! This routine performs a single analysis step for
! PDAF in offline mode using PDAF with domain-decomposition.
!
! The analysis is performed by calling a filter-specific 
! routine PDAF\_put\_state\_X.
!
! In this routine, the real names of most of the 
! user-supplied routines for PDAF are specified (see below).
!
! !REVISION HISTORY:
! 2009-11 - Lars Nerger - Initial code by restructuring
! Later revisions - see svn log
!
! !USES:
  USE mod_parallel_pdaf, &    ! Parallelization
       ONLY: mype_world, abort_parallel
  USE mod_assimilation_pdaf, & ! Variables for assimilation
       ONLY: filtertype

  IMPLICIT NONE

! ! External subroutines 
! !  (subroutine names are passed over to PDAF in the calls to 
! !  PDAF_get_state and PDAF_assimilate_X. This allows the user 
! !  to specify the actual name of a routine. However, the 
! !  PDAF-internal name of a subroutine might be different from
! !  the external name!)
!
  ! Interface between model and PDAF, and prepoststep
  EXTERNAL :: collect_state_pdaf, &    ! Routine to collect a state vector from model fields
       prepoststep_ens_offline         ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, &  ! Provide number of local analysis domains
       init_dim_l_pdaf, &             ! Initialize state dimension for local analysis domain
       g2l_state_pdaf, &              ! Get state on local analysis domain from global state
       l2g_state_pdaf                 ! Update global state from state on local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: &
       init_dim_obs_pdafomi, &         ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &               ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi          ! Get dimension of obs. vector for local analysis domain

! !CALLING SEQUENCE:
! Called by: main
! Calls: PDAF_get_state (possible, but not required!)
! Calls: PDAF_put_state_seik
! Calls: PDAF_put_state_enkf
! Calls: PDAF_put_state_lseik
! Calls: PDAF_put_state_etkf
! Calls: PDAF_put_state_letkf
! Calls: PDAF_put_state_lenkf
! Calls: PDAF_put_state_netf
! Calls: PDAF_put_state_lnetf
! Calls: MPI_barrier (MPI)
!EOP

! local variables
  INTEGER :: status    ! Status flag for filter routines


! ************************
! *** Perform analysis ***
! ************************

! *** Note on PDAF_get_state for offline implementation: ***
! *** For the offline mode of PDAF the call to           ***
! *** PDAF_get_state is not required as no forecasting   ***
! *** is performed in this mode. However, it is save     ***
! *** to call PDAF_get_state, even it is not necessary.  ***
! *** The functionality of PDAF_get_state is deactived   ***
! *** for the offline mode.                              ***


  IF (filtertype == 7) THEN
     CALL PDAFomi_put_state_local(collect_state_pdaf, init_dim_obs_pdafomi, &
          obs_op_pdafomi, prepoststep_ens_offline, init_n_domains_pdaf, init_dim_l_pdaf, &
          init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, status)
  ELSEIF (filtertype == 5) THEN
     WRITE(*,*) 'NOTICE: NO ANALYSIS STEP IS EXECUTED FOR filtertype=5'
     status = 0
  ELSE
     WRITE(*,*) 'ERROR in init_pdaf: Only LESKTF (filtertype=7) is implemented'
     status = -1
  END IF


! ************************
! *** Check error flag ***
! ************************

  IF (status /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a47,i4,a1/)') &
          'ERROR ', status, &
          ' during assimilation with PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

END SUBROUTINE assimilation_pdaf_offline
