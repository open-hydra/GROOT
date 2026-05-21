module GROOT_Register_Discretization
  implicit none
  private
  public :: Register_Discretization

contains

  subroutine Register_Discretization()
    use GROOT_Input_Registry
    use GROOT_Config_Types_m, only: obj_dtm
    use GROOT_Global_m
    implicit none
    character(len=:), allocatable :: sec
    sec = trim(codename)//'-Discretization'

    call reg%add(sec, 'Nr', obj_dtm%Nr, '256', &
        'Number of ray directions (uniform sphere sampling)', '> 0', .true.)

    call reg%add(sec, 'tol', obj_dtm%tol, '1.0e-3', &
        'Convergence tolerance on the wall-flux residual', '> 0', .true.)

    call reg%add(sec, 'itermax', obj_dtm%itermax, '10', &
        'Maximum number of DTM iterations', '>= 0', .false.)

    call reg%add(sec, 'ang', obj_dtm%wedge_angle_deg, '1.0', &
        'Wedge half-angle [deg] for 2-D axisymmetric cases (twoDax = .true.)', '> 0', .false.)

  end subroutine Register_Discretization

end module GROOT_Register_Discretization
