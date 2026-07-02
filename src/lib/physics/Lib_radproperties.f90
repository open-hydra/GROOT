module GROOT_Lib_Radproperties
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  integer(kind=I4) :: iH2O, iCO2, iCO

  contains

  ! finds index of radiative species in cfd species array; loads SNBW database when needed
  subroutine map_species
    use GROOT_Global_m, only: model, cfd_species, nsc
    use GROOT_Lib_SNBW, only: read_snb_param
    use GROOT_SNB_Datadir, only: snb_data_path
    implicit none

    if (model == 'wsgg-H2O') then
      iH2O = findindex('H2O')
    elseif (model == 'wsgg-H2OCO2') then
      iH2O = findindex('H2O')
      iCO2 = findindex('CO2')
    elseif (model == 'snbw' .or. model == 'snb') then
      iH2O = findindex('H2O')
      iCO2 = findindex('CO2')
      iCO  = findindex('CO')
      call read_snb_param(snb_data_path)
    elseif (model == 'slw') then
      iH2O = findindex('H2O')
      iCO2 = findindex('CO2')
      iCO  = findindex('CO')
    else

    endif

  end subroutine map_species

  ! finds which entry of cfd species array corresponds to species "spec"
  function findindex(spec) result(ind)
    use GROOT_Global_m, only: cfd_species, nsc
    implicit none 
    character(len=*), intent(in) :: spec 
    integer :: ind, is 
    ind = 0

    do is = 1, nsc
      if (trim(cfd_species(is))==trim(spec)) then 
        ind = is 
        exit 
      endif
    enddo 
  end function findindex

  ! computes radiative properties of mixture
  subroutine compute_gasproperties(xvec,p,T,fvs,ka,a,Ib,kb)
    use GROOT_Global_m, only: model, nsc, nscdat, ind_species, cfd_species, sigma, pi, k_user, Ngg
    use GROOT_Lib_SNBW, only: co_ksnb, co_ksnb_malkmus
    use GROOT_Lib_SLW,  only: co_kslw
    implicit none
    real(kind=R8), intent(in)           :: xvec(1:nsc)  !< molar fractions
    real(kind=R8), intent(in)           :: p            !< pressure [Pa]
    real(kind=R8), intent(in)           :: T            !< temperature [K]
    real(kind=R8), intent(in)           :: fvs          !< soot volume fraction (reserved)
    real(kind=R8), intent(out)          :: ka(0:Ngg)    !< absorption coefficient [m-1]
    real(kind=R8), intent(out)          :: a(0:Ngg)     !< spectral weight [-]
    real(kind=R8), intent(out)          :: Ib           !< total blackbody intensity [W m-2 sr-1]
    real(kind=R8), optional, intent(out):: kb(0:)       !< Malkmus beta_eff (SNB full only)

    integer(kind=I4) :: is
    real(kind=R8) :: eps, xH2O, xCO2, xCO_loc
    real(kind=R8) :: beta_loc(0:Ngg)

    if (model=='const') then
      ka = k_user
      a  = 1.0_R8
      Ib = sigma*T**4/pi
    elseif (model == 'gray') then
      ka = 0.0_R8
      do is = 1, nscdat
        if (ind_species(is).ne.0) then
          ka = ka + xvec(ind_species(is))*graygas(T,p,is)
        endif
      enddo
      Ib = sigma*T**4/pi
      a  = 1.0_R8
    elseif (model == 'wsgg-H2O') then
      xH2O = 0.0_R8
      if (iH2O > 0) xH2O = xvec(iH2O)
      call wsggH2O(xH2O, p, T, 1.0_R8, ka, a, eps)
      Ib = sigma*T**4/pi
    elseif (model == 'wsgg-H2OCO2') then
      xH2O = 0.0_R8
      xCO2 = 0.0_R8
      if (iH2O > 0) xH2O = xvec(iH2O)
      if (iCO2 > 0) xCO2 = xvec(iCO2)
      call wsggH2OCO2_MR(xH2O, xCO2, p, T, 1.0_R8, ka, a, eps)
      Ib = sigma*T**4/pi
    elseif (model == 'snbw') then
      xH2O    = 0.0_R8
      xCO2    = 0.0_R8
      xCO_loc = 0.0_R8
      if (iH2O > 0) xH2O    = xvec(iH2O)
      if (iCO2 > 0) xCO2    = xvec(iCO2)
      if (iCO  > 0) xCO_loc = xvec(iCO)
      call co_ksnb(T, p, xH2O, xCO2, xCO_loc, ka, a)
      Ib = sigma*T**4/pi
    elseif (model == 'snb') then
      xH2O    = 0.0_R8
      xCO2    = 0.0_R8
      xCO_loc = 0.0_R8
      if (iH2O > 0) xH2O    = xvec(iH2O)
      if (iCO2 > 0) xCO2    = xvec(iCO2)
      if (iCO  > 0) xCO_loc = xvec(iCO)
      call co_ksnb_malkmus(T, p, xH2O, xCO2, xCO_loc, ka, a, beta_loc)
      if (present(kb)) kb = beta_loc
      Ib = sigma*T**4/pi
    elseif (model == 'slw') then
      xH2O    = 0.0_R8
      xCO2    = 0.0_R8
      xCO_loc = 0.0_R8
      if (iH2O > 0) xH2O    = xvec(iH2O)
      if (iCO2 > 0) xCO2    = xvec(iCO2)
      if (iCO  > 0) xCO_loc = xvec(iCO)
      call co_kslw(T, p, xH2O, xCO2, xCO_loc, Ngg, ka, a)
      Ib = sigma*T**4/pi
    else
      print *, 'selected model not available: ', trim(model)
      stop
    endif
  end subroutine compute_gasproperties

!> @brief Function computes the absorption coefficient
!> H2O, CO2 and CO data are taken from R.S. Barlow, A.N. Karpetis, J.H. Frank, J.-Y. Chen, Scalar profiles and no formation
!in laminar opposed-flow partially premixed methane/air flames, Combust. Flame
!127 (3) (2001) 2102–2118, https://doi.org/10.1016/S0010-2180(01)00313-3 for T< 2500 K, with an hyperbolic extrapolation 
!for T>2500 fitted on the data reported in rivere 2012 https://doi.org/10.1016/j.ijheatmasstransfer.2012.03.019
!> CH4 data from chu 2014 https://doi.org/10.1007/S11708-013-0292-4
!> C2H4 data from tuntomo 1989 http://dx.doi.org/10.1080/08916158908946356
!> RP1 data is taken as the arithmetic average of the k values given in NISTIR 6646, considering a refractive index of 1.44
!> HCL data is taken from S. P. Fuss, Determination of Planck Mean Absorption Coefficients for HBr,HCL, and HF
!> Journal of Heat Transfer 2002, https://doi.org/10.1115/1.1416689
  function graygas(T,p,is) result(kval)
    implicit none
    real(kind=R8),    intent(in) :: T    !< Temperature in K
    real(kind=R8),    intent(in) :: p    !< Pressure in Pa
    integer(kind=I4), intent(in) :: is   !< Species index as in sc_database
    real(kind=R8)             :: kval  ! absorption coefficient

    real(kind=R8) :: T0, T1, KP0, DKP0, a(1:5)
    
    select case (is)
    case(1) !H20 
      if (T<2500.) then
        kval = -0.23093+1E+03*(-1.12390)*(T**(-1))+1E+06*(9.41530)*(T**(-2))+1E+09*(-2.99880)*(T**(-3))+1E+12*(0.51382)*(T**(-4))&
        +1E+15*(-1.86840E-05)*(T**(-5)) 
      else
        T0=2500.
        T1=0.
        KP0 = -0.23093+1E+03*(-1.12390)*(T0**(-1))+1E+06*(9.41530)*(T0**(-2))+1E+09*(-2.99880)*(T0**(-3))+1E+12*(0.51382)*(T0**(-4))&
        +1E+15*(-1.86840E-05)*(T**(-5)) 
        DKP0 = -1E+03*(-1.12390)*(T0**(-2))-2E+06*( 9.41530)*(T0**(-3)) -3E+09*(-2.99880)*(T0**(-4)) -4E+12*(0.51382)*(T0**(-5))&
        -5E+15*(-1.86840E-05)*(T0**(-6))
        kval = KP0*(((T-T1)/(T0-T1))**(((DKP0/KP0)*(T0-T1))))
      endif
      kval = kval*0.9869*p/1.d5
    case(2) !CO2
      if (T<2500.) then
        kval = 18.741 +1E+03*(-123.310)*(T**(-1)) +1E+06*( 273.500)*(T**(-2)) +1E+09*(-194.050)*(T**(-3)) &
        +1E+12*(  56.310)*(T**(-4))+1E+15*( -5.8169)*(T**(-5))
      else
        T0= 2500.
        T1=-2500.
        KP0=18.741 +1E+03*(-123.310)*(T0**(-1)) +1E+06*( 273.500)*(T0**(-2)) +1E+09*(-194.050)*(T0**(-3)) &
        +1E+12*(  56.310)*(T0**(-4))+1E+15*( -5.8169)*(T0**(-5))
        DKP0 = -1E+03*(-123.310)*(T0**(-2)) -2E+06*( 273.500)*(T0**(-3)) -3E+09*(-194.050)*(T0**(-4)) &
        -4E+12*(  56.310)*(T0**(-5)) -5E+15*( -5.8169)*(T0**(-6))
        kval = KP0*(((T-T1)/(T0-T1))**((DKP0/KP0)*(T0-T1)))
      endif
      kval = kval*0.9869*p/1.d5
    case default
      ! error message
    end select
 
  end function graygas



!> @brief Subroutine containing WSGG model by Fabiani et al. (JPP 2025) for H2O
subroutine wsggH2O(xH2O, p, TC, L, kp, a, eps)
  implicit none
  real(kind=R8), intent(in) :: xH2O         !< Molar fractions
  real(kind=R8), intent(in) :: p           !< Pressure
  real(kind=R8), intent(in) :: TC           !< Temperature
  real(kind=R8), intent(in) :: L           !< Path lenght
  real(kind=R8), intent(out):: kp(0:4)
  real(kind=R8), intent(out):: a(0:4)
  real(kind=R8), intent(out):: eps

  real(kind=R8) :: cij(1:4,1:4),k(0:4), Tref
  integer(kind=I4) :: ig, j 

  Tref = 2400.d0 

  k(0) = 0.d0
  k(1:4) = (/  0.967803,   6.310879,   0.156390,   0.014594/) 
  cij(1,1:4) = (/ -0.080801,   1.163064,  -1.160097,   0.314726/) 
  cij(2,1:4) = (/  1.009603,  -1.917783,   1.244422,  -0.274227/) 
  cij(3,1:4) = (/ -0.306869,   1.640970,  -1.278695,   0.292618/) 
  cij(4,1:4) = (/  0.120045,   0.126245,   0.027891,  -0.040793/) 

  eps = 0.d0
  do ig=1,4
    a(ig) = 0
    do j=1,4
      a(ig) = a(ig)+cij(ig,j)*(TC/Tref)**(j-1)
    enddo
    kp(ig) = k(ig)*xH2O*p/1.d5
    eps=eps+a(ig)*(1-exp(-kp(ig)*L))
  enddo
  a(0)=1-sum(a(1:4))
  
end subroutine wsggH2O



!> @brief Subroutine containing WSGG model by Fabiani et al. (JPP 2025) for H2O-CO2 mixtures
subroutine wsggH2OCO2_MR(xH2O, xCO2, pc, TC, L, kp, a, eps)
  implicit none
  real(kind=R8), intent(in) :: pc           !< Pressure
  real(kind=R8), intent(in) :: TC           !< Temperature
  real(kind=R8), intent(in) :: xH2O,xCO2    !< Molar fractions
  real(kind=R8), intent(in) :: L            !< Path lenght
  real(kind=R8), intent(out):: kp(0:4)
  real(kind=R8), intent(out):: a(0:4)
  real(kind=R8), intent(out):: eps

  real(kind=R8) :: MR,k(1:4),cij(1:4,1:4), Tref, MRvec(1:11),MRarr(1:11),T
  integer(kind=I4) :: iMR, ig, j

  Tref = 2300.d0
  if (TC<1000.d0) then
    T = 1000.d0
  else 
    T = TC 
  endif

  MR = xH2O/(xCO2+1.d-5)
  MRvec = (/0.125d0,0.25d0,0.5d0,0.75d0,1.d0,2.d0,2.5d0,3d0,4.d0,6.d0,8.d0/)
  MRarr = abs(MRvec-MR)
  iMR = minloc(MRarr,1) 

  select case (iMR)
  case(1)
    k(1:4) = (/  0.201355,  83.681211,   2.177752,   0.018921/) 
    cij(1,1:4) = (/  0.338670,   0.126308,  -0.373209,   0.123367/) 
    cij(2,1:4) = (/  0.198925,  -0.301228,   0.165138,  -0.031874/) 
    cij(3,1:4) = (/  0.582898,  -0.935571,   0.565397,  -0.120346/) 
    cij(4,1:4) = (/ -0.532591,   2.300377,  -1.761992,   0.402114/) 
    case(2)
    k(1:4) = (/  0.176943,  68.107015,   1.927536,   0.019695/) 
    cij(1,1:4) = (/  0.013091,   0.948018,  -0.939157,   0.244841/) 
    cij(2,1:4) = (/  0.252767,  -0.430343,   0.267167,  -0.058113/) 
    cij(3,1:4) = (/  0.728348,  -1.080919,   0.581898,  -0.108590/) 
    cij(4,1:4) = (/ -0.396877,   1.790107,  -1.318230,   0.291051/) 
    case(3)
    k(1:4) = (/  0.167276,   1.952047,  52.441380,   0.018473/) 
    cij(1,1:4) = (/ -0.206884,   1.486154,  -1.253509,   0.298121/) 
    cij(2,1:4) = (/  0.725780,  -0.884108,   0.345353,  -0.036283/) 
    cij(3,1:4) = (/  0.335381,  -0.627902,   0.419616,  -0.096346/) 
    cij(4,1:4) = (/ -0.236388,   1.262742,  -0.904787,   0.195389/) 
    case(4)
    k(1:4) = (/  0.169308,  42.244993,   1.969425,   0.017736/) 
    cij(1,1:4) = (/ -0.277568,   1.632683,  -1.308007,   0.298863/) 
    cij(2,1:4) = (/  0.401490,  -0.778331,   0.531477,  -0.123635/) 
    cij(3,1:4) = (/  0.670017,  -0.660894,   0.137143,   0.020843/) 
    cij(4,1:4) = (/ -0.159840,   1.029515,  -0.728327,   0.155358/) 
    case(5)
    k(1:4) = (/  0.174844,  36.680098,   1.995313,   0.017628/) 
    cij(1,1:4) = (/ -0.302002,   1.667088,  -1.300099,   0.289894/) 
    cij(2,1:4) = (/  0.450271,  -0.887988,   0.611855,  -0.142996/) 
    cij(3,1:4) = (/  0.617410,  -0.486243,  -0.015639,   0.061352/) 
    cij(4,1:4) = (/ -0.120589,   0.914228,  -0.640886,   0.135256/) 
    case(6)
    k(1:4) = (/  0.186016,  25.467792,   1.950971,   0.017696/) 
    cij(1,1:4) = (/ -0.320818,   1.637206,  -1.207716,   0.254019/) 
    cij(2,1:4) = (/  0.584577,  -1.174975,   0.814065,  -0.190224/) 
    cij(3,1:4) = (/  0.458295,  -0.017773,  -0.401982,   0.160087/) 
    cij(4,1:4) = (/ -0.056075,   0.724775,  -0.497380,   0.102112/) 
    case(7)
    k(1:4) = (/  0.189185,  23.112324,   1.928168,   0.017814/) 
    cij(1,1:4) = (/ -0.318111,   1.608209,  -1.169021,   0.241602/) 
    cij(2,1:4) = (/  0.625255,  -1.259805,   0.872198,  -0.203453/) 
    cij(3,1:4) = (/  0.406785,   0.124317,  -0.514596,   0.188062/) 
    cij(4,1:4) = (/ -0.043536,   0.687196,  -0.468620,   0.095378/) 
    case(8)
    k(1:4) = (/  0.191435,  21.511787,   1.908193,   0.017914/) 
    cij(1,1:4) = (/ -0.314168,   1.581880,  -1.137337,   0.231845/) 
    cij(2,1:4) = (/  0.656536,  -1.324673,   0.916222,  -0.213366/) 
    cij(3,1:4) = (/  0.366303,   0.234028,  -0.600378,   0.209141/) 
    cij(4,1:4) = (/ -0.035419,   0.662399,  -0.449544,   0.090890/) 
    case(9)
    k(1:4) = (/  0.194536,  19.480684,   1.878634,   0.018078/) 
    cij(1,1:4) = (/ -0.305586,   1.538113,  -1.089016,   0.217523/) 
    cij(2,1:4) = (/  0.701315,  -1.417862,   0.979208,  -0.227454/) 
    cij(3,1:4) = (/  0.306826,   0.392967,  -0.723220,   0.239031/) 
    cij(4,1:4) = (/ -0.026015,   0.632872,  -0.426646,   0.085466/) 
    case(10)
    k(1:4) = (/  0.198014,  17.336022,   1.842568,   0.018294/) 
    cij(1,1:4) = (/ -0.290310,   1.475621,  -1.025942,   0.199634/) 
    cij(2,1:4) = (/  0.755432,  -1.533041,   1.057727,  -0.245063/) 
    cij(3,1:4) = (/  0.232303,   0.590153,  -0.874267,   0.275481/) 
    cij(4,1:4) = (/ -0.018249,   0.606954,  -0.406312,   0.080618/) 
    case(11)
    k(1:4) = (/  0.199856,  16.157754,   1.820514,   0.018416/) 
    cij(1,1:4) = (/ -0.278339,   1.432731,  -0.985239,   0.188458/) 
    cij(2,1:4) = (/  0.788192,  -1.605694,   1.108504,  -0.256651/) 
    cij(3,1:4) = (/  0.185639,   0.713984,  -0.969253,   0.298410/) 
    cij(4,1:4) = (/ -0.015392,   0.596299,  -0.397852,   0.078594/) 
    end select

    eps = 0.d0
    do ig=1,4
      a(ig) = 0
      do j=1,4
        a(ig) = a(ig)+cij(ig,j)*(T/Tref)**(j-1)
      enddo
      kp(ig) = k(ig)*(xH2O+xCO2)*pc/1.d5
      eps=eps+a(ig)*(1-exp(-kp(ig)*L))
    enddo
    a(0)=1-sum(a(1:4))
 
end subroutine wsggH2OCO2_MR


end module GROOT_Lib_Radproperties