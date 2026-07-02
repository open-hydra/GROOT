!> @brief RC-SLW spectral model via the radlib C++ library (BYUignite/radlib).
!>
!> Rank-Correlated Spectral Line-based Weighted-sum-of-gray-gases (RC-SLW).
!> radlib returns nGGa = nGG+1 values (nGG gray gases + 1 clear gas), indexed 0..nGG:
!>   kabs(0) = 0        -> transparent window (WSGG-compatible convention)
!>   awts(0..nGG)       -> weights, sum to 1
!> so GROOT's Ngg (extra gray gases beyond the transparent index 0) maps directly to
!> radlib's nGG, and the (0:Ngg) arrays line up with no resizing.
!>
!> The rad_rcslw object fixes P and TbTref at construction; TbTref must equal the local
!> gas temperature passed to get_k_a, so we construct one object per call.
!> Species: H2O, CO2, CO (CH4 accepted but ignored internally by radlib).
!>
!> This file is compiled only when USE_RADLIB=ON; otherwise Lib_slw_stub.f90 provides
!> a same-signature stub that aborts with an explanatory message.
module GROOT_Lib_SLW
  use iso_c_binding
  implicit none
  private
  public :: co_kslw

  ! radlib RC-SLW tabulated pressure range [atm]; outside it radlib calls exit(0),
  ! which would abort the whole solver. We clamp instead and warn once.
  real(kind=8), parameter :: SLW_PMIN_PA = 0.1_8  * 101325.0_8   ! 0.1 atm
  real(kind=8), parameter :: SLW_PMAX_PA = 50.0_8 * 101325.0_8   ! 50  atm
  logical, save :: p_clamp_warned = .false.

  ! --- radlib C interface (src/fortran/c_interface.cc); scalars by reference (no VALUE) ---
  interface
    function rad_rcslw_C_interface(nGG, TbTref, P, fvsoot, xH2O, xCO2, xCO) &
             result(rad_ptr) bind(C, name="rad_rcslw_C_interface")
      import :: c_ptr, c_int, c_double
      type(c_ptr)    :: rad_ptr
      integer(c_int) :: nGG
      real(c_double) :: TbTref, P, fvsoot, xH2O, xCO2, xCO
    end function rad_rcslw_C_interface

    subroutine get_k_a_C_interface(rad_ptr, kabs, awts, T, P, fvsoot, &
                                   xH2O, xCO2, xCO, xCH4) &
               bind(C, name="get_k_a_C_interface")
      import :: c_ptr, c_double
      type(c_ptr), value           :: rad_ptr
      real(c_double), dimension(*) :: kabs, awts
      real(c_double)               :: T, P, fvsoot, xH2O, xCO2, xCO, xCH4
    end subroutine get_k_a_C_interface

    subroutine rad_delete_C_interface(rad_ptr) bind(C, name="rad_delete_C_interface")
      import :: c_ptr
      type(c_ptr), value :: rad_ptr
    end subroutine rad_delete_C_interface
  end interface

contains

  !> @brief Fill ka(0:Ngg_slw) and a(0:Ngg_slw) for the RC-SLW model.
  !> Same output layout as co_ksnb / wsgg routines (index 0 = transparent, ka=0).
  subroutine co_kslw(T, p, xH2O, xCO2, xCO, Ngg_slw, ka, a)
    real(kind=8), intent(in)  :: T, p          !< temperature [K], pressure [Pa]
    real(kind=8), intent(in)  :: xH2O, xCO2, xCO  !< molar fractions [-]
    integer,      intent(in)  :: Ngg_slw       !< number of gray gases (= radlib nGG)
    real(kind=8), intent(out) :: ka(0:Ngg_slw) !< absorption coefficient [m-1]
    real(kind=8), intent(out) :: a (0:Ngg_slw) !< spectral weight [-]

    type(c_ptr)    :: obj
    real(c_double) :: kabs_c(0:Ngg_slw), awts_c(0:Ngg_slw)  ! nGG+1 slots (radlib writes 0..nGG)
    real(c_double) :: Tc, Pc, xH2Oc, xCO2c, xCOc
    real(c_double) :: zero
    real(kind=8)   :: p_use
    integer(c_int) :: ngg_c

    ! Clamp pressure to radlib's tabulated range so it never calls exit(0).
    p_use = min(max(p, SLW_PMIN_PA), SLW_PMAX_PA)
    if (p_use /= p .and. .not. p_clamp_warned) then
      write(*,'(A,ES10.3,A)') ' WARNING (SLW): pressure out of radlib range [0.1,50] atm, '// &
        'clamping to ', p_use/101325.0_8, ' atm (further warnings suppressed).'
      p_clamp_warned = .true.
    end if

    zero   = 0.0_c_double
    ngg_c  = int(Ngg_slw, c_int)
    Tc     = real(T,     c_double)
    Pc     = real(p_use, c_double)   ! radlib constructor converts Pa -> atm internally
    xH2Oc  = real(xH2O,  c_double)
    xCO2c  = real(xCO2, c_double)
    xCOc   = real(xCO,  c_double)

    ! TbTref = Tc (defines the ALBDF F-grid boundaries); must match the T passed to get_k_a.
    obj = rad_rcslw_C_interface(ngg_c, Tc, Pc, zero, xH2Oc, xCO2c, xCOc)
    call get_k_a_C_interface(obj, kabs_c, awts_c, Tc, Pc, zero, xH2Oc, xCO2c, xCOc, zero)
    call rad_delete_C_interface(obj)

    ka(0:Ngg_slw) = real(kabs_c(0:Ngg_slw), kind=8)
    a (0:Ngg_slw) = real(awts_c(0:Ngg_slw), kind=8)
    ka(0) = 0.0_8   ! transparent window (radlib already sets kabs(0)=0; enforce explicitly)
  end subroutine co_kslw

end module GROOT_Lib_SLW
