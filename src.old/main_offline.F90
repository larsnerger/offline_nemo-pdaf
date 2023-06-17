!$Id$
!BOP
!
! !ROUTINE: main --- Main program for example of PDAF offline implementation
!
! !INTERFACE:
program MAIN_OFFLINE

! !DESCRIPTION:
! This is the main program for an example implementation of
! PDAF with domain-decomposition and offline configuration.
!
! In the offline mode, we assume that the ensemble
! integrations are performed in a separate program (model)
! and the forecasted ensemble can be read from files. After
! initializing the ensemble information by reading model
! outputs, a single analysis step is performed. Subsequently,
! the analysis ensemble can be written to files that can be 
! used to initialize another ensemble forecast.
!
! Using PDAF for domain-decomposition, the offline
! mode can be used to perform assimilation with domain-
! decomposed models. If the models write results for each 
! sub-domain, these can be read here using the same 
! parallelization. Then, the filter analysis can be 
! performed utilizing this parallelization. If the files
! contain the full model state, PDAF in offline mode
! can be used either on a single processor, or the 
! fields can be distributed in this program to utilize
! the parallelization of the filters.
!
! Parameters can be set in the code, or - preferably -
! by command line arguments that are parsed by the 
! module PARSER. The format for this is
! EXECUTABLE -HANDLE1 VALUE1 -HANDLE2 VALUE2 ...
! The handles are defined in the code before the calls
! to the routine PARSE.
!
! !REVISION HISTORY:
! 2008-07 - Lars Nerger - Initial code
! Later revisions - see svn log
!
! !USES:
  use mod_parallel_pdaf, &     ! Parallelization
       only: MPI_COMM_WORLD, MPIerr, npes_world, mype_world, &
       init_parallel, finalize_parallel
  use mod_assimilation_pdaf, & ! Variables for assimilation
       only: program_mode
  use mod_memcount_pdaf, &
       only: memcount_ini, memcount_get
  use timer, &
       only: timeit, time_tot

  implicit none
!EOP


! Local variables
  integer :: i                 ! Counter


! **********************
! *** Initialize MPI ***
! **********************

  call init_parallel() ! initializes MPI


! ********************************
! ***      INITIALIZATION      ***
! ********************************

  ! Initialize memory counting and timers
  call memcount_ini(4)
  call timeit(4, 'ini')

  call timeit(1,'new')


! *** Initial Screen output ***
  initscreen: if (mype_world == 0) then

     write (*, '(/8x, a/)') '+++++ PDAF offline mode +++++'
     write (*, '(9x, a)') 'Data assimilation with PDAF'

     if (npes_world > 1) then
        write (*, '(/16x, a, i3, a/)') 'Running on ', npes_world, ' PEs'
     else
        write (*, '(/16x, a/)') 'Running on 1 PE'
     end if
     write (*, '(/)')
     
  end if initscreen

  
! *** Initialize MPI communicators for PDAF (model and filter) ***
! *** NOTE: It is always n_modeltasks=1 for offline mode       ***

  call init_parallel_pdaf(0, 1)

! *** Initialize model information ***
! *** This should only be information on the model dimension
! *** Generally, this could be joined with init_pdaf.

  call timeit(2,'new')
  call initialize()
  call timeit(2,'old')


! *******************************
! ***      ASSIMILATION       ***
! *******************************

  ! *** Initialize PDAF ***

  call timeit(3,'new')
  call init_pdaf()
  call timeit(3,'old')


  ! *** Perform analysis ***

  if (trim(program_mode)=='assim') then
     if (mype_world == 0) &
          write (*, '(/2x, a)') 'PDAF offline mode: START ASSIMILATION'

     call timeit(4,'new')
     call assimilation_pdaf_offline()
     call timeit(4,'old')

  elseif (trim(program_mode)=='covar') then
     if (mype_world == 0) &
          write (*, '(/2x, a)') 'PDAF offline: GENERATE COVARIANCE MATRIX'

     call timeit(4,'new')
     call eofcovar()
     call timeit(4,'old')

  end if


  ! Synchronize at barrier for exit
  call MPI_Barrier(MPI_COMM_WORLD, MPIerr) 
  write (*,*) 'model PE exited: mype ', mype_world

  call timeit(1,'old')


! ********************
! *** Finishing up ***
! ********************

! *** Final screen output ***
  if (mype_world == 0) then
       write (*, '(/1x, a)') 'PDAF offline mode: EXITED ASSIMILATION'

     ! *** Show allocated memory for the model ***
     write (*, '(/22x, a)') 'Model - Memory overview'
     write (*, '(14x, 45a)') ('-', i=1, 45)
     write (*, '(25x, a)') 'Allocated memory  (MB)'
     write (*, '(18x, a, f12.3, a)') &
          'Model field:', memcount_get(1, 'M'), ' MB (persistent)'

     write (*, '(16x, a, f12.3, a)') &
          'ensemble init:', memcount_get(2, 'M'), ' MB (temporary)'
     write (*, '(17x, a, f12.3, a)') &
          'Pre-Poststep:', memcount_get(3, 'M'), ' MB (temporary)'
     write (*, '(20x, a, f12.3, a)') &
          'gen_covar:', memcount_get(4, 'M'), ' MB (temporary)'
  end if

  ! *** Finalize PDAF - print memory and timing information
  call finalize_pdaf(0)

  if (mype_world == 0) then
     write (*, '(/17x, a)') 'Offline - Timing information'
     write (*, '(10x, 45a)') ('-', i=1, 45)
     ! Timing summary for assimilation
     write (*, '(19x, a, F11.3, 1x, a)') 'initialize model:', time_tot(2), 's'
     write (*, '(18x, a, F11.3, 1x, a)') 'initialize filter:', time_tot(3), 's'
     if (trim(program_mode)=='assim') then
        write (*, '(23x, a, F11.3, 1x, a)') 'assimilation:', time_tot(4), 's'
     else
        write (*, '(27x, a, F11.3, 1x, a)') 'gen_covar:', time_tot(4), 's'
     end if
     write (*, '(19x, a, F11.3, 1x, a)') 'total run time:', time_tot(1), 's'

     write (*, '(/1x, a)') 'PDAF offline mode: END'
  end if


! *** Terminate MPI
  call finalize_parallel()

end program MAIN_OFFLINE
