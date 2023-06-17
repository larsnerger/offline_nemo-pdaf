!$Id$
!BOP
!
! !ROUTINE: init_pdaf_info - Screen output on assimilation configuration
!
! !INTERFACE:
SUBROUTINE init_pdaf_info()

! !DESCRIPTION:
! This routine performs a model-sided screen output about
! the coniguration of the data assimilation system.
! Using this output is optional. Most of the information
! is also displayed by PDAF itself when it is initialized
! in PDAF_init. Not displayed by PDAF is the assimilation
! interval (delt_obs), which is unknown to PDAF.
!
! !REVISION HISTORY:
! 2011-05 - Lars Nerger - Initial code extracted from init_pdaf
! Later revisions - see svn log
!
! !USES:
  USE mod_assimilation_pdaf, & ! Variables for assimilation
       ONLY: screen, filtertype, subtype, dim_ens, delt_obs, &
       model_error, model_err_amp, incremental, covartype, &
       type_forget, forget, epsilon, rank_analysis_enkf, locweight, &
       int_rediag, type_trans

  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: init_pdaf
!EOP


! *****************************
! *** Initial Screen output ***
! *****************************

  IF (filtertype == 0) THEN
     WRITE (*, '(a,/21x, a)') 'NEMO-PDAF','Filter: SEEK'
     IF (subtype == 2) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- fixed basis filter with update of matrix U'
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- no re-diagonalization of VUV^T'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- fixed basis filter & no update of matrix U'
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- no re-diagonalization of VUV^T'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Offline mode'
     END IF
     WRITE (*, '(a,13x, a, i5)') 'NEMO-PDAF','number of EOFs:', dim_ens
     IF (subtype /= 5) WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF','Assimilation interval:', delt_obs
     WRITE (*, '(a,10x, a, f5.2)') 'NEMO-PDAF','forgetting factor:', forget
     IF (subtype /= 5) THEN
        IF ((int_rediag > 0) .AND. ((subtype /= 2) .OR. (subtype /= 3))) &
             WRITE (*, '(a,10x, a, i4, a)') 'NEMO-PDAF',&
             'Re-diag each ', int_rediag, '-th analysis step'
     ELSE
        IF (int_rediag == 1) THEN
           WRITE (*, '(a,10x, a)') 'NEMO-PDAF','Perform re-diagonalization'
        ELSE
           WRITE (*, '(a,10x, a)') 'NEMO-PDAF','No re-diagonalization'
        END IF
     END IF
  ELSE IF (filtertype == 1) THEN
     WRITE (*, '(a,21x, a)') 'NEMO-PDAF','Filter: SEIK'
     IF (subtype == 2) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- fixed state covariance matrix'
     ELSE IF (subtype == 4) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- use ensemble transformation'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Offline mode'
     END IF
     WRITE (*, '(a,14x, a, i5)') 'NEMO-PDAF','ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF','Assimilation interval:', delt_obs
     WRITE (*, '(a,10x, a, f5.2)') 'NEMO-PDAF','forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(a,6x, a, f5.2)') 'NEMO-PDAF','model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 2) THEN
     WRITE (*, '(a,21x, a)') 'NEMO-PDAF','Filter: EnKF'
     IF (subtype == 5) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Offline mode'
     END IF
     WRITE (*, '(a,14x, a, i5)') 'NEMO-PDAF','ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF','Assimilation interval:', delt_obs
     WRITE (*, '(a,10x, a, f5.2)') 'NEMO-PDAF','forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(a,6x, a, f5.2)') 'NEMO-PDAF','model error amplitude:', model_err_amp
     END IF
     IF (rank_analysis_enkf > 0) THEN
        WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF',&
             'analysis with pseudo-inverse of HPH, rank:', rank_analysis_enkf
     END IF
  ELSE IF (filtertype == 3) THEN
     WRITE (*, '(a,21x, a)') 'NEMO-PDAF','Filter: LSEIK'
     IF (subtype == 2) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- fixed error-space basis'
     ELSE IF (subtype == 3) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- fixed state covariance matrix'
     ELSE IF (subtype == 4) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- use ensemble transformation'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Offline mode'
     END IF
     WRITE (*, '(a,14x, a, i5)') 'NEMO-PDAF','ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF','Assimilation interval:', delt_obs
     WRITE (*, '(a,10x, a, f5.2)') 'NEMO-PDAF','forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(a,6x, a, f5.2)') 'NEMO-PDAF','model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 4) THEN
     WRITE (*, '(a,21x, a)') 'NEMO-PDAF','Filter: ETKF'
     IF (subtype == 0) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Variant using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Variant following Hunt et al. (2007)'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Offline mode'
     END IF
     WRITE (*, '(a,14x, a, i5)') 'NEMO-PDAF','ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF','Assimilation interval:', delt_obs
     WRITE (*, '(a,10x, a, f5.2)') 'NEMO-PDAF','forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(a,6x, a, f5.2)') 'NEMO-PDAF','model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 5) THEN
     WRITE (*, '(a,21x, a)') 'NEMO-PDAF','Filter: LETKF'
     IF (subtype == 0) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Variant using T-matrix'
     ELSE IF (subtype == 1) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Variant following Hunt et al. (2007)'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Offline mode'
     END IF
     WRITE (*, '(a,14x, a, i5)') 'NEMO-PDAF','ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF','Assimilation interval:', delt_obs
     WRITE (*, '(a,10x, a, f5.2)') 'NEMO-PDAF','forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(a,6x, a, f5.2)') 'NEMO-PDAF','model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 6) THEN
     WRITE (*, '(a,21x, a)') 'NEMO-PDAF','Filter: ESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Offline mode'
     END IF
     WRITE (*, '(a,14x, a, i5)') 'NEMO-PDAF','ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF','Assimilation interval:', delt_obs
     WRITE (*, '(a,10x, a, f5.2)') 'NEMO-PDAF','forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*,'(a,6x, a, f5.2)') 'NEMO-PDAF','model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 7) THEN
     WRITE (*, '(a,21x, a)') 'NEMO-PDAF','Filter: LESTKF'
     IF (subtype == 0) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Standard mode'
     ELSE IF (subtype == 5) THEN
        WRITE (*, '(a,6x, a)') 'NEMO-PDAF','-- Offline mode'
     END IF
     WRITE (*, '(a,14x, a, i5)') 'NEMO-PDAF','ensemble size:', dim_ens
     IF (subtype /= 5) WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF','Assimilation interval:', delt_obs
     WRITE (*, '(a,10x, a, f5.2)') 'NEMO-PDAF','forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(a,6x, a, f5.2)') 'NEMO-PDAF','model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 8) THEN
     WRITE (*, '(a,21x, a)') 'NEMO-PDAF','Filter: LEnKF'
     WRITE (*, '(a,14x, a, i5)') 'NEMO-PDAF','ensemble size:', dim_ens
     WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF','Assimilation interval:', delt_obs
     WRITE (*, '(a,10x, a, f5.2)') 'NEMO-PDAF','forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(a,6x, a, f5.2)') 'NEMO-PDAF','model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 9) THEN
     WRITE (*, '(a,21x, a)') 'NEMO-PDAF','Filter: NETF'
     WRITE (*, '(a,14x, a, i5)') 'NEMO-PDAF','ensemble size:', dim_ens
     WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF','Assimilation interval:', delt_obs
     WRITE (*, '(a,10x, a, f5.2)') 'NEMO-PDAF','forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(a,6x, a, f5.2)') 'NEMO-PDAF','model error amplitude:', model_err_amp
     END IF
  ELSE IF (filtertype == 10) THEN
     WRITE (*, '(a,21x, a)') 'NEMO-PDAF','Filter: LNETF'
     WRITE (*, '(a,14x, a, i5)') 'NEMO-PDAF','ensemble size:', dim_ens
     WRITE (*, '(a,6x, a, i5)') 'NEMO-PDAF','Assimilation interval:', delt_obs
     WRITE (*, '(a,10x, a, f5.2)') 'NEMO-PDAF','forgetting factor:', forget
     IF (model_error) THEN
        WRITE (*, '(a,6x, a, f5.2)') 'NEMO-PDAF','model error amplitude:', model_err_amp
     END IF
  END IF     

END SUBROUTINE init_pdaf_info
