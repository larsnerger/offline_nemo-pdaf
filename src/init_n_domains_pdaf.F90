!$Id$
!> Set number of local analysis domains
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF
!!
!! The routine is called in PDAF_X_update 
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains. 
!! It has to set the number of local analysis 
!! domains for the PE-local domain.
!!
!! Implementation for the 2D online example
!! without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
subroutine init_n_domains_pdaf(step, n_domains_p)

   use mod_nemo_pdaf, &
        only: nwet

  implicit none

! *** Arguments ***
  integer, intent(in)  :: step        !< Current time step
  integer, intent(out) :: n_domains_p !< PE-local number of analysis domains


! ************************************
! *** Initialize number of domains ***
! ************************************
 
   n_domains_p = nwet

end subroutine init_n_domains_pdaf
