module GROOT_Register_IO
  use iso_fortran_env, only: R8 => real64
  implicit none
  private
  public :: Register_IO


contains

  subroutine Register_IO()
    use GROOT_Input_Registry
    use GROOT_Config_Types_m, only: obj_io
    use GROOT_Global_m
    implicit none
    character(len=:), allocatable :: sec

    sec = 'MOSE-IO'

    call reg%add( sec, 'sol-format', obj_io%gas_format, 'tecplot ascii', 'Gas Solution (OUTPUT/field.*) format', &
     'tecplot ascii, tecplot binary, vtk ascii, vtk raw', .false. )

    sec = 'GROOT-IO'

    call reg%add( sec, 'sol-format', obj_io%sol_format, 'tecplot ascii', 'Radiative Solution (OUTPUT/rad-field.*) format', &
     'tecplot ascii, tecplot binary, vtk ascii, vtk raw', .false. )   

  end subroutine Register_IO



end module GROOT_Register_IO
