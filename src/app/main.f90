program GROOT_program
#if defined (_OPENMP)
  use omp_lib
#endif
  use GROOT_Advanced_Types_m, only: GROOT_simulation_type
  use GROOT_Config_Types_m,   only: obj_sim_param
  use GROOT_Procedures_m,     only: GROOT_type
  implicit none
  type(GROOT_type)           :: GROOT
  type(GROOT_simulation_type):: simulation

#if defined (_OPENMP)
  !$omp parallel
  obj_sim_param%nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  write(*,'(A)')    ' Parallel execution'
  write(*,'(A)')    ' OpenMP:'
  write(*,'(A,I4)') ' Number of threads --> ', obj_sim_param%nthreads
#else
  write(*,'(A)')    ' Serial execution'
  obj_sim_param%nthreads = 1
#endif

  call GROOT%setup(simulation)

  call GROOT%solve(simulation)

  call GROOT%postprocess(simulation)

contains 

  subroutine Dummy_Function
    ! Empty subroutine to be passed as an argument to GROOT%solve. 
    ! It can be used for user-defined operations during the solution process.
    ! It is used in HYDRA coupling procedures.
  end subroutine Dummy_Function


end program GROOT_program
