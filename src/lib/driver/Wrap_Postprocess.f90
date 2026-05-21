module GROOT_Wrap_Postprocess
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: GROOT_postprocess

  character(69), private :: shell_format = "('GROOT | Time of operation:', F10.3, ' min')"

contains

  subroutine GROOT_postprocess(simulation)
    use GROOT_Advanced_Types_m, only: GROOT_simulation_type
    use GROOT_Config_Types_m,   only: obj_sim_param, obj_io
    use GROOT_Global_m,         only: GROOT_phase_prefix
    use GROOT_IO_Solution,      only: write_solution
    use GROOT_IO_Wall,          only: Initialize_Wall_File, Write_Wall_Solution
    implicit none
    type(GROOT_simulation_type), intent(inout) :: simulation
    real(R8) :: sim_time

    !! Timing
    call cpu_time(obj_sim_param%cputime(2))
    sim_time = ( obj_sim_param%cputime(2) - obj_sim_param%cputime(1) ) / obj_sim_param%nthreads

    write(*,*)
    write(*, shell_format) sim_time / 60.0_R8

    !! Write radiative field solution
    call write_solution(simulation%domain, simulation%ORP, &
                        trim(GROOT_phase_prefix)//'field')

    !! Write wall solution
    call Initialize_Wall_File(simulation%domain)
    call Write_Wall_Solution(simulation%domain, trim(GROOT_phase_prefix)//'wall', obj_io%sol_format)

  end subroutine GROOT_postprocess

end module GROOT_Wrap_Postprocess
