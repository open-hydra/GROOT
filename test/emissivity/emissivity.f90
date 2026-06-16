!! ============================================================
!! Homogeneous column emissivity test.
!!
!! Computes total emissivity of a uniform H2O-CO2 mixture at
!! fixed composition and path length over a grid of (T, p) values
!! using three spectral models:
!!   - WSGG   : Fabiani et al. JPP 2025, 4 gray gases
!!   - SNBW   : SNB weak limit (Beer-Lambert per band)
!!   - SNB    : SNB full Malkmus transmissivity
!!
!! Fixed parameters:
!!   xH2O = 0.2,  xCO2 = 0.1,  xCO = 0.0
!!   L    = 1.0 m
!!
!! Output: OUTPUT/emissivity.dat
!!   Columns: T[K]  p[bar]  eps_wsgg  eps_snbw  eps_snb
!! ============================================================
program emissivity_test
  use iso_fortran_env,          only: R8 => real64
  use GROOT_Global_m,           only: sigma, pi
  use GROOT_Lib_SNBW,           only: N_SNB_BANDS, read_snb_param, &
                                      co_ksnb, co_ksnb_malkmus
  use GROOT_SNB_Datadir,        only: snb_data_path
  use GROOT_Lib_Radproperties,  only: wsggH2OCO2_MR
  implicit none

  ! ---- fixed mixture and geometry ----
  real(R8), parameter :: xH2O = 0.2_R8
  real(R8), parameter :: xCO2 = 0.1_R8
  real(R8), parameter :: xCO  = 0.0_R8
  real(R8), parameter :: L    = 1.0_R8   ! path length [m]

  ! ---- parameter grids ----
  integer,  parameter :: NT = 6, NP = 5
  real(R8) :: T_vals(NT) = [1000.0_R8, 1500.0_R8, 2000.0_R8, &
                             2500.0_R8, 3000.0_R8, 3500.0_R8]  ! [K]
  real(R8) :: p_bar (NP) = [1.0_R8, 10.0_R8, 50.0_R8, 100.0_R8, 300.0_R8] ! [bar]

  ! ---- local work arrays ----
  real(R8) :: ka  (0:N_SNB_BANDS)
  real(R8) :: a   (0:N_SNB_BANDS)
  real(R8) :: beta(0:N_SNB_BANDS)
  real(R8) :: kp_wsgg(0:4), a_wsgg(0:4)

  real(R8) :: T, p_pa, tau
  real(R8) :: eps_wsgg, eps_snbw, eps_snb
  integer  :: iT, ip, ig, unit_out

  ! ---- load SNB database ----
  call read_snb_param(snb_data_path)

  ! ---- prepare output directory ----
  call execute_command_line('mkdir -p OUTPUT', wait=.true.)

  open(newunit=unit_out, file='OUTPUT/emissivity.dat', action='write')
  write(unit_out,'(A)') '# Homogeneous column emissivity: xH2O=0.2  xCO2=0.1  xCO=0.0  L=1m'
  write(unit_out,'(A)') '#'//repeat('-',65)
  write(unit_out,'(A)') '#   T[K]   p[bar]    eps_wsgg    eps_snbw    eps_snb'
  write(unit_out,'(A)') '#'//repeat('-',65)

  do iT = 1, NT
    T = T_vals(iT)
    do ip = 1, NP
      p_pa = p_bar(ip) * 1.0e5_R8   ! bar → Pa

      ! ---- WSGG ----
      call wsggH2OCO2_MR(xH2O, xCO2, p_pa, T, L, kp_wsgg, a_wsgg, eps_wsgg)

      ! ---- SNBW : Beer-Lambert per band ----
      call co_ksnb(T, p_pa, xH2O, xCO2, xCO, ka, a)
      eps_snbw = 0.0_R8
      do ig = 1, N_SNB_BANDS
        eps_snbw = eps_snbw + a(ig) * (1.0_R8 - exp(-ka(ig) * L))
      end do

      ! ---- SNB full : Malkmus transmissivity per band ----
      call co_ksnb_malkmus(T, p_pa, xH2O, xCO2, xCO, ka, a, beta)
      eps_snb = 0.0_R8
      do ig = 1, N_SNB_BANDS
        if (ka(ig) > 1.0e-30_R8) then
          tau = exp(-beta(ig) * (sqrt(1.0_R8 + 2.0_R8*ka(ig)*L/beta(ig)) - 1.0_R8))
        else
          tau = 1.0_R8
        end if
        eps_snb = eps_snb + a(ig) * (1.0_R8 - tau)
      end do

      write(unit_out,'(F8.1, F8.1, 3F12.6)') T, p_bar(ip), eps_wsgg, eps_snbw, eps_snb
    end do
    write(unit_out,*)
  end do

  close(unit_out)

  write(*,'(A)') ''
  write(*,'(A)') ' Emissivity table written to OUTPUT/emissivity.dat'
  write(*,'(A)') ' Columns: T[K]  p[bar]  eps_wsgg  eps_snbw  eps_snb'
  write(*,'(A)') ''

end program emissivity_test
