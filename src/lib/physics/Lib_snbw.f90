!> @brief Statistical Narrow Band model — weak limit (SNBW) and full Malkmus (SNB).
!>
!> Each of the N_SNB_BANDS = 450 narrow bands (25 cm-1 wide, covering 25..11250 cm-1).
!> SNBW treats each band as Beer-Lambert gray gas: tau = exp(-k_bar * L).
!> SNB uses the Malkmus transmissivity: tau = exp(-beta*(sqrt(1+2*k_bar*u/beta)-1))
!> where u = x*P_atm*L_cm is the column amount and beta is the Malkmus parameter.
!>
!> Database: Soufiani & Taine (1997) for H2O and CO2; Riviere & Soufiani (2012) for CO.
!> Species included: H2O, CO2, CO.  CH4 and HCL are not supported.
!> xO2 is set to zero in the H2O gamma broadening formula (standard for H2O-CO2-N2 mixtures).

module GROOT_Lib_SNBW
  use, intrinsic :: iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private

  public :: N_SNB_BANDS
  public :: read_snb_param
  public :: co_ksnb
  public :: co_ksnb_malkmus
  public :: snb_planck_weight

  !> Number of 25 cm-1 narrow bands (eta = 25, 50, ..., 11250 cm-1)
  integer, parameter :: N_SNB_BANDS = 450

  ! SNB database arrays: K = mean abs coeff [cm-1 atm-1], D = 1/delta [cm] (inverse line spacing)
  !   temperature grid: 300, 400, ..., 5000 K  (48 points, step 100 K)
  !   band grids are species-specific (see FINDI)
  real(R8) :: KH2O(48, 449), DH2O(48, 449)
  real(R8) :: KCO2(48, 323), DCO2(48, 323)
  real(R8) :: KCO (48, 194), DCO (48, 194)

contains

  ! -----------------------------------------------------------------------
  !> @brief Band-integrated Planck weight for narrow band ig at temperature T.
  !>
  !> Computes  a(ig) = pi * B_ig(T) * Delta_eta_m / (sigma * T^4)
  !> where B_ig [W m-2 sr-1 m] is the spectral blackbody function at
  !> the band centre eta = 25*ig [cm-1] and Delta_eta_m = 2500 m-1.
  !> The factor pi is already included in the planck() function.
  !> Summing over ig=1..N_SNB_BANDS recovers (approximately) sigma*T^4/pi,
  !> i.e.  sum(a) ≈ 1.
  elemental function snb_planck_weight(T, ig) result(a_ig)
    use GROOT_Global_m, only: sigma
    real(R8), intent(in) :: T   !< temperature [K]
    integer,  intent(in) :: ig  !< band index 1..N_SNB_BANDS
    real(R8) :: a_ig
    real(R8), parameter :: DELTA_ETA_M = 2500.0_R8  ! 25 cm-1 converted to m-1
    a_ig = planck(T, 25.0_R8 * ig, 1.0_R8) * DELTA_ETA_M / (sigma * T**4)
  end function snb_planck_weight

  ! -----------------------------------------------------------------------
  !> @brief Per-band absorption coefficient [m-1] and Planck weights for H2O + CO2 + CO.
  !>
  !> Fills ka(0:N_SNB_BANDS) and a(0:N_SNB_BANDS).
  !> Band 0 is the WSGG-style transparent window; for SNBW it is set to zero
  !> (no transparent window — all spectral energy is carried by bands 1..450).
  subroutine co_ksnb(T, P, xH2O, xCO2, xCO, ka, a)
    implicit none
    real(R8), intent(in)  :: T                        !< temperature [K]
    real(R8), intent(in)  :: P                        !< pressure [Pa]
    real(R8), intent(in)  :: xH2O, xCO2, xCO         !< molar fractions [-]
    real(R8), intent(out) :: ka(0:N_SNB_BANDS)        !< absorption coefficient [m-1]
    real(R8), intent(out) :: a (0:N_SNB_BANDS)        !< Planck fraction [-]

    integer  :: iwvnb, ih2o, ico2, ico, iT
    real(R8) :: wvnb, RT
    logical  :: LIH2O, LICO2, LICO

    ! band 0: transparent window placeholder (unused in SNBW, kept for interface compatibility)
    ka(0) = 0.0_R8
    a (0) = 0.0_R8

    call TMNO(T, RT, iT)   ! T is constant across bands; compute once

    wvnb = 0.0_R8
    do iwvnb = 1, N_SNB_BANDS
      wvnb = 25.0_R8 * iwvnb   ! band centre [cm-1]

      call FINDI(wvnb, ico2, ih2o, ico, LICO2, LIH2O, LICO)

      ka(iwvnb) = 0.0_R8
      if (LICO2) ka(iwvnb) = ka(iwvnb) + xCO2 * (KCO2(iT,iCO2) + RT*(KCO2(iT+1,iCO2) - KCO2(iT,iCO2)))
      if (LIH2O) ka(iwvnb) = ka(iwvnb) + xH2O * (KH2O(iT,iH2O) + RT*(KH2O(iT+1,iH2O) - KH2O(iT,iH2O)))
      if (LICO)  ka(iwvnb) = ka(iwvnb) + xCO  * (KCO (iT,iCO)  + RT*(KCO (iT+1,iCO)  - KCO (iT,iCO) ))

      ! Database units: [cm-1 atm-1].  Convert to [m-1]:
      !   × (P/1e5) converts Pa → bar (~atm),  × 100 converts cm-1 → m-1
      ka(iwvnb) = ka(iwvnb) * (P / 1.0d5) * 100.0_R8

      ! Planck weight: fraction of total blackbody power in this band
      a(iwvnb) = snb_planck_weight(T, iwvnb)
    end do

  end subroutine co_ksnb

  ! -----------------------------------------------------------------------
  !> @brief Per-band ka [m-1], Planck weights, and Malkmus beta_eff for H2O + CO2 + CO.
  !>
  !> Fills ka(0:N_SNB_BANDS), a(0:N_SNB_BANDS), beta_eff(0:N_SNB_BANDS).
  !> beta_eff(ig) is the dimensionless Malkmus parameter for band ig, defined as
  !>   beta_eff = sum_sp(x_sp * P_atm * k_sp * beta_sp) / sum_sp(x_sp * P_atm * k_sp)
  !> where beta_sp = 2 * gamma_sp * D_sp  (gamma [cm-1], D = 1/delta [cm]).
  !> This is independent of path length and can be precomputed per cell.
  !> The transmissivity along a ray segment of length ds [m] is then:
  !>   tau(ig) = exp(-beta_eff * (sqrt(1 + 2*ka(ig)*ds/beta_eff) - 1))
  !> Broadening formulas from TRSMI (Riviere & Soufiani 2012); xO2 set to zero.
  subroutine co_ksnb_malkmus(T, P, xH2O, xCO2, xCO, ka, a, beta_eff)
    implicit none
    real(R8), intent(in)  :: T                        !< temperature [K]
    real(R8), intent(in)  :: P                        !< pressure [Pa]
    real(R8), intent(in)  :: xH2O, xCO2, xCO         !< molar fractions [-]
    real(R8), intent(out) :: ka      (0:N_SNB_BANDS)  !< absorption coefficient [m-1]
    real(R8), intent(out) :: a       (0:N_SNB_BANDS)  !< Planck fraction [-]
    real(R8), intent(out) :: beta_eff(0:N_SNB_BANDS)  !< Malkmus beta parameter [-]

    integer  :: iwvnb, ih2o, ico2, ico, iT
    real(R8) :: wvnb, RT
    real(R8) :: P_atm, T296
    real(R8) :: k_CO, D_CO, GAM_CO, beta_CO
    real(R8) :: k_CO2, D_CO2, GAM_CO2, beta_CO2
    real(R8) :: k_H2O, D_H2O, GAM_H2O, beta_H2O
    real(R8) :: SK, SB, w
    logical  :: LIH2O, LICO2, LICO

    ka(0)       = 0.0_R8
    a(0)        = 0.0_R8
    beta_eff(0) = 1.0e10_R8

    P_atm = P / 1.0d5    ! Pa → bar (consistent with co_ksnb and WSGG)
    T296  = 296.0_R8 / T

    call TMNO(T, RT, iT)   ! T is constant across bands; compute once

    do iwvnb = 1, N_SNB_BANDS
      wvnb = 25.0_R8 * iwvnb
      call FINDI(wvnb, ico2, ih2o, ico, LICO2, LIH2O, LICO)

      SK = 0.0_R8
      SB = 0.0_R8

      if (LICO) then
        k_CO   = KCO(iT,ico)  + RT*(KCO(iT+1,ico)  - KCO(iT,ico) )
        D_CO   = DCO(iT,ico)  + RT*(DCO(iT+1,ico)  - DCO(iT,ico) )
        GAM_CO = P_atm * (0.075_R8*xCO2*T296**0.6_R8 &
                        + 0.06_R8*(1.0_R8-xCO2-xH2O)*T296**0.7_R8 &
                        + 0.12_R8*xH2O*T296**0.82_R8)
        beta_CO = 2.0_R8 * GAM_CO * D_CO
        w  = xCO * P_atm * k_CO
        SK = SK + w;  SB = SB + w * beta_CO
      end if

      if (LICO2) then
        k_CO2   = KCO2(iT,ico2) + RT*(KCO2(iT+1,ico2) - KCO2(iT,ico2))
        D_CO2   = DCO2(iT,ico2) + RT*(DCO2(iT+1,ico2) - DCO2(iT,ico2))
        GAM_CO2 = P_atm * (0.07_R8*xCO2 + 0.058_R8*(1.0_R8-xCO2-xH2O) &
                         + 0.10_R8*xH2O) * T296**0.7_R8
        beta_CO2 = 2.0_R8 * GAM_CO2 * D_CO2
        w   = xCO2 * P_atm * k_CO2
        SK  = SK + w;  SB = SB + w * beta_CO2
      end if

      if (LIH2O) then
        k_H2O   = KH2O(iT,ih2o) + RT*(KH2O(iT+1,ih2o) - KH2O(iT,ih2o))
        D_H2O   = DH2O(iT,ih2o) + RT*(DH2O(iT+1,ih2o) - DH2O(iT,ih2o))
        ! xO2 = 0 (not tracked; standard for H2O-CO2-N2 combustion mixtures)
        GAM_H2O = P_atm * (0.462_R8*T296*xH2O &
                         + sqrt(T296)*(0.0792_R8*(1.0_R8-xCO2) + 0.106_R8*xCO2))
        beta_H2O = 2.0_R8 * GAM_H2O * D_H2O
        w    = xH2O * P_atm * k_H2O
        SK   = SK + w;  SB = SB + w * beta_H2O
      end if

      ka(iwvnb) = SK * 100.0_R8   ! [cm-1] → [m-1]
      a(iwvnb)  = snb_planck_weight(T, iwvnb)
      if (SK > 1.0e-30_R8) then
        beta_eff(iwvnb) = SB / SK
      else
        beta_eff(iwvnb) = 1.0e10_R8  ! transparent band → Beer-Lambert limit
      end if
    end do

  end subroutine co_ksnb_malkmus

  ! -----------------------------------------------------------------------
  !> @brief Read SNB spectral database from datadir.
  !>
  !> File format (Soufiani & Taine / Riviere & Soufiani):
  !>   SNBCO  — 194 K records, then 194 D records
  !>   SNBCO2 — 323 K records, then 323 D records
  !>   SNBH2O — 449 K records, then 449 D records
  !>   Each record: wavenumber  val_1  val_2 ... val_48   (48 temperature points)
  !>   K [cm-1 atm-1], D [cm] = 1/(mean line spacing delta [cm-1])
  subroutine read_snb_param(datadir)
    implicit none
    character(len=*), intent(in) :: datadir

    integer  :: i, j
    real(R8) :: xxnu

    open(unit=45, file=trim(datadir)//'SNBCO',  status='old', action='read')
    open(unit=46, file=trim(datadir)//'SNBCO2', status='old', action='read')
    open(unit=47, file=trim(datadir)//'SNBH2O', status='old', action='read')

    ! --- CO (194 bands) ---
    do i = 1, 194
      read(45,*) xxnu, (KCO(j,i), j=1,48)
    end do
    do i = 1, 194
      read(45,*) xxnu, (DCO(j,i), j=1,48)
    end do
    close(45)

    ! --- CO2 (323 bands) ---
    do i = 1, 323
      read(46,*) xxnu, (KCO2(j,i), j=1,48)
    end do
    do i = 1, 323
      read(46,*) xxnu, (DCO2(j,i), j=1,48)
    end do
    close(46)

    ! --- H2O (449 bands) ---
    do i = 1, 449
      read(47,*) xxnu, (KH2O(j,i), j=1,48)
    end do
    do i = 1, 449
      read(47,*) xxnu, (DH2O(j,i), j=1,48)
    end do
    close(47)

  end subroutine read_snb_param

  ! -----------------------------------------------------------------------
  !> @brief Spectral band index for each species at wavenumber wvnb [cm-1].
  subroutine FINDI(wvnb, ico2, ih2o, ico, LICO2, LIH2O, LICO)
    implicit none
    real(R8), intent(in)  :: wvnb
    integer,  intent(out) :: ico2, ih2o, ico
    logical,  intent(out) :: LICO2, LIH2O, LICO

    ! CO2: 250..8300 cm-1, 323 bands of 25 cm-1
    if (wvnb >= 250.0_R8 .and. wvnb <= 8300.0_R8) then
      ico2 = int((wvnb - 250.0_R8) / 25.0_R8) + 1
    else
      ico2 = -9
    end if

    ! H2O: 50..11250 cm-1, 449 bands of 25 cm-1
    if (wvnb >= 50.0_R8 .and. wvnb <= 11250.0_R8) then
      ih2o = int((wvnb - 50.0_R8) / 25.0_R8) + 1
    else
      ih2o = -9
    end if

    ! CO: 1600..6425 cm-1, 194 bands of 25 cm-1
    if (wvnb >= 1600.0_R8 .and. wvnb <= 6425.0_R8) then
      ico = int((wvnb - 1600.0_R8) / 25.0_R8) + 1
    else
      ico = -9
    end if

    LICO2 = (ico2 > 0)
    LIH2O = (ih2o > 0)
    LICO  = (ico  > 0)
  end subroutine FINDI

  ! -----------------------------------------------------------------------
  !> @brief Temperature interpolation index and fractional offset.
  !>
  !> SNB temperature grid: 300, 400, ..., 5000 K  (48 points, step 100 K).
  !> Returns iT in [1,47] and RT in [0,1) such that
  !>   k(T) ≈ K(iT,band) + RT*(K(iT+1,band) - K(iT,band))
  subroutine TMNO(TM, RT, iT)
    implicit none
    real(R8), intent(in)  :: TM
    real(R8), intent(out) :: RT
    integer,  intent(out) :: iT

    if (TM > 300.0_R8) then
      if (TM < 5000.0_R8) then
        RT = (TM - 300.0_R8) / 100.0_R8
        iT = int(RT + 1.0e-6_R8)
        RT = RT - iT
        iT = iT + 1
      else
        RT = 1.0_R8
        iT = 47
      end if
    else
      RT = 0.0_R8
      iT = 1
    end if
  end subroutine TMNO

  ! -----------------------------------------------------------------------
  !> @brief Spectral blackbody function (already multiplied by pi) [W m-2 m].
  !>
  !> Returns  pi * B_eta(T, n)  in [W m-2 (m-1)-1]
  !> so that  integral_{0}^{inf} planck(T,eta,n) d(eta_m) = sigma*T^4 / n^2.
  elemental function planck(T, eta, n) result(E)
    real(R8), intent(in) :: T    !< temperature [K]
    real(R8), intent(in) :: eta  !< wavenumber [cm-1]
    real(R8), intent(in) :: n    !< refractive index
    real(R8) :: E
    real(R8), parameter :: C1 = 3.7418d-16  ! W m2  (= 2hc^2)
    real(R8), parameter :: C2 = 1.4388d0    ! cm K  (= hc/k in cm·K)
    E = C1 * (eta * 1.0d2)**3 / (n**2 * (exp(C2 * (eta + 1.0d-6) / (n * T)) - 1.0_R8))
  end function planck

end module GROOT_Lib_SNBW
