!> @brief Stub for the RC-SLW model when radlib is not compiled in (USE_RADLIB=OFF).
!>
!> Provides the same public interface as the real GROOT_Lib_SLW module so the rest of
!> GROOT compiles unchanged.  Any attempt to actually run model=slw aborts with a
!> message telling the user to reconfigure with -DUSE_RADLIB=ON.
module GROOT_Lib_SLW
  implicit none
  private
  public :: co_kslw

contains

  subroutine co_kslw(T, p, xH2O, xCO2, xCO, Ngg_slw, ka, a)
    real(kind=8), intent(in)  :: T, p
    real(kind=8), intent(in)  :: xH2O, xCO2, xCO
    integer,      intent(in)  :: Ngg_slw
    real(kind=8), intent(out) :: ka(0:Ngg_slw), a(0:Ngg_slw)

    ! silence unused-argument warnings, then abort
    ka = 0.0_8
    a  = 0.0_8
    write(*,'(A)') ' ERROR: model=slw requires radlib. Reconfigure with cmake -DUSE_RADLIB=ON.'
    stop 1
  end subroutine co_kslw

end module GROOT_Lib_SLW
