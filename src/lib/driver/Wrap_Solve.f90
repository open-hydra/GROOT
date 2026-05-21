module GROOT_Wrap_Solve
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: GROOT_solve

contains

  subroutine GROOT_solve(simulation)
    use GROOT_Advanced_Types_m, only: GROOT_simulation_type
    use GROOT_Config_Types_m,   only: obj_dtm
    use GROOT_Lib_DTM,          only: co_flux, co_source
    implicit none
    type(GROOT_simulation_type), intent(inout) :: simulation

    !! Compute radiative heat flux on all boundary faces
    call co_flux(simulation%domain)

    !! Optionally compute volumetric source term (div q)
    if (obj_dtm%source) call co_source(simulation%domain, 2 * obj_dtm%Nr)

  end subroutine GROOT_solve

end module GROOT_Wrap_Solve
