module GROOT_Register_Options
  implicit none
  private
  public :: Register_Options

contains

  subroutine Register_Options()
    use GROOT_Input_Registry
    use GROOT_Config_Types_m, only: obj_dtm
    use GROOT_Global_m
    implicit none
    character(len=:), allocatable :: sec
    sec = trim(codename)//'-Options'

    call reg%add(sec, 'source', obj_dtm%source, '.false.', &
        'Compute volumetric radiative source term div(q_r)', '.true. , .false.', .false.)

    call reg%add(sec, 'twoDax', obj_dtm%twoDax, '.false.', &
        '2-D axisymmetric mode: domain is a thin wedge in the x-r plane', '.true. , .false.', .false.)

    call reg%add(sec, 'optically_thin', obj_dtm%optically_thin, '.false.', &
        'Optically-thin approximation for source term', '.true. , .false.', .false.)

    
  end subroutine Register_Options

end module GROOT_Register_Options
