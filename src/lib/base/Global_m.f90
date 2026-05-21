module GROOT_Global_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use GROOT_Parameters_m   ! re-exports hlen, llen, clen, codenames, prefixes

  implicit none

  !! ----------------------------------------------------------
  !! Physical constants
  real(R8) :: pi    = 4._8 * datan(1._8)
  real(R8) :: sigma = 5.6703d-8
  real(R8), parameter :: Ru_ = 8314.51d0

  !! Runtime globals — set by Assign_Setup after reading ini
  integer :: Nr                   !< Number of rays
  integer :: Ngg                  !< Number of gray gases (or narrow bands)
  integer :: itermax              !< Max DTM iterations
  integer :: imax = 10000         !< Max ray-marching steps
  real(R8) :: wedge_angle_deg      !< Wedge angle for 2-D axisymmetric (deg)
  real(R8) :: tol                  !< DTM heat-flux convergence tolerance


  !! ----------------------------------------------------------
  !! Radiation model flags — set by Assign_Setup
  real(R8)         :: eps_wall
  real(R8)         :: k_user
  character(llen)  :: model
  logical          :: source
  logical          :: kcost
  logical          :: twoDax
  logical          :: optically_thin

  !! ----------------------------------------------------------
  !! Species parameters — set by Read_IdealGas_Properties
  integer                              :: nsc, nscrad, nscdat
  character(len=clen), allocatable     :: cfd_species(:)
  character(len=clen), allocatable     :: rad_species(:)
  integer,             allocatable     :: ind_species(:)
  real(R8),            allocatable     :: wm_tab(:)

  !! ----------------------------------------------------------
  !! Auxiliary pointers into wall.tec variable list
  integer :: nT    = 1
  integer :: nsoot = 0
  integer :: nrans = 0

end module GROOT_Global_m
