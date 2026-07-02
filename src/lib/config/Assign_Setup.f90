module GROOT_Assign_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Assign_Setup

contains

  subroutine Assign_Setup()
    use GROOT_Config_Types_m, only: obj_dtm, obj_rad_model
    use GROOT_Global_m
    implicit none

    !! DTM options
    Nr             = obj_dtm%Nr
    tol            = obj_dtm%tol
    itermax        = obj_dtm%itermax
    wedge_angle_deg= obj_dtm%wedge_angle_deg
    source         = obj_dtm%source
    twoDax         = obj_dtm%twoDax
    optically_thin = obj_dtm%optically_thin

    !! Radiation model
    model          = obj_rad_model%model
    eps_wall       = obj_rad_model%eps_wall
    k_user         = obj_rad_model%k_user
    p_ref          = obj_rad_model%p_ref

    !! Constant absorption coefficient flag
    kcost = (trim(model) == 'const')

    !! Number of gray gases / narrow bands
    select case (trim(model))
      case ('gray', 'const');            Ngg = 0
      case ('wsgg-H2O', 'wsgg-H2OCO2'); Ngg = 4
      case ('snbw', 'snb');              Ngg = 450  ! N_SNB_BANDS in Lib_snbw
      case ('slw')
        N_SLW_GASES = obj_rad_model%slw_ngray       ! radlib nGG (gray gases, clear gas excluded)
        Ngg         = N_SLW_GASES
      case default;                      Ngg = 0
    end select

  end subroutine Assign_Setup

end module GROOT_Assign_Setup
