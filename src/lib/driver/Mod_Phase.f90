module GROOT_Mod_Phase
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private
  public :: Setup_Gas
  public :: Setup_Rad
  public :: Setup_Wall
  public :: Setup_Face
  public :: Setup_SNBWall
  public :: Setup_WSGGWall
contains

  subroutine Setup_Gas(grid, IOfield_gas)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use Lib_ORION_Data
    use GROOT_IO_Solution, only: IOtime
    use GROOT_Global_m, only: nsc, nrans, nsoot
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    type(orion_data),        intent(inout) :: IOfield_gas
    integer :: b


    do b = 1, grid%nb
      associate(Nx => grid%blk(b)%dim(1), &
                Ny => grid%blk(b)%dim(2), &
                Nz => grid%blk(b)%dim(3))
        grid%blk(b)%gas_phase%rho(1:nsc, 1:Nx, 1:Ny, 1:Nz) = &
          IOfield_gas%block(b)%vars(1:nsc, :, :, :)
        grid%blk(b)%gas_phase%p(1:Nx, 1:Ny, 1:Nz) = &
          IOfield_gas%block(b)%vars(nsc+4, :, :, :)
        grid%blk(b)%gas_phase%T(1:Nx, 1:Ny, 1:Nz) = &
          IOfield_gas%block(b)%vars(nsc+nrans+nsoot+5, :, :, :)
        grid%blk(b)%gas_phase%R(1:Nx, 1:Ny, 1:Nz) = &
          IOfield_gas%block(b)%vars(nsc+nrans+nsoot+7, :, :, :)
      end associate
    end do

    grid%time = IOtime
  end subroutine Setup_Gas


  subroutine Setup_Rad(grid, IOfield_rad)
    use GROOT_Lib_Radproperties, only: compute_gasproperties, map_species
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use Lib_ORION_Data
    use GROOT_Global_m, only: nsc, Ngg, model
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    type(orion_data),        intent(out)   :: IOfield_rad
    integer :: i, j, k, b
    real(R8) :: ka(0:Ngg), ks_loc(0:Ngg), a(0:Ngg), Ib, fvs, kb(0:Ngg)
    logical  :: snb_full

    call Compute_Molfractions(grid)
    call map_species()

    snb_full = (trim(model) == 'snb')

    allocate(IOfield_rad%block(1:grid%nb))
    do b = 1, grid%nb
      associate(Nx => grid%blk(b)%dim(1), &
                Ny => grid%blk(b)%dim(2), &
                Nz => grid%blk(b)%dim(3))
        IOfield_rad%block(b)%Ni = Nx
        IOfield_rad%block(b)%Nj = Ny
        IOfield_rad%block(b)%Nk = Nz

        allocate(IOfield_rad%block(b)%vars(1:4, 1:Nx, 1:Ny, 1:Nz))
        allocate(IOfield_rad%block(b)%mesh(1,   1:Nx, 1:Ny, 1:Nz))
        allocate(grid%blk(b)%ka    (1:Nx, 1:Ny, 1:Nz, 0:Ngg))
        allocate(grid%blk(b)%a     (1:Nx, 1:Ny, 1:Nz, 0:Ngg))
        allocate(grid%blk(b)%Ib    (1:Nx, 1:Ny, 1:Nz))
        allocate(grid%blk(b)%source(1:Nx, 1:Ny, 1:Nz))
        allocate(grid%blk(b)%G     (1:Nx, 1:Ny, 1:Nz))
        if (snb_full) allocate(grid%blk(b)%kb(1:Nx, 1:Ny, 1:Nz, 0:Ngg))

        do k = 1, Nz; do j = 1, Ny; do i = 1, Nx
          fvs = 0.0_R8

          if (snb_full) then
            call compute_gasproperties( &
              xvec = grid%blk(b)%gas_phase%rho(1:nsc, i, j, k), &
              p    = grid%blk(b)%gas_phase%p(i, j, k),           &
              T    = grid%blk(b)%gas_phase%T(i, j, k),           &
              fvs  = fvs, ka = ka, a = a, Ib = Ib, kb = kb)
            grid%blk(b)%kb(i,j,k, 0:Ngg) = kb
          else
            call compute_gasproperties( &
              xvec = grid%blk(b)%gas_phase%rho(1:nsc, i, j, k), &
              p    = grid%blk(b)%gas_phase%p(i, j, k),           &
              T    = grid%blk(b)%gas_phase%T(i, j, k),           &
              fvs  = fvs, ka = ka, a = a, Ib = Ib)
          end if
          if (isnan(sum(ka)) .or. isnan(sum(a))) then
            write(*, '(A,4I4)') ' NaN in radiative coefficients at block/cell:', b, i, j, k
            stop
          end if
          grid%blk(b)%ka    (i,j,k, 0:Ngg) = ka
          grid%blk(b)%a     (i,j,k, 0:Ngg) = a
          grid%blk(b)%Ib    (i,j,k)        = Ib
          grid%blk(b)%source(i,j,k)        = 0.0_R8
          grid%blk(b)%G     (i,j,k)        = 0.0_R8
        end do; end do; end do
      end associate
    end do
  end subroutine Setup_Rad


  !! @brief Set up wall temperatures and emissivities from wall.tec file.
  subroutine Setup_Wall(grid, IOfield_wall)
    use strings,       only: parse
    use GROOT_IO_Wall, only: viscous_flag
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use Lib_ORION_Data
    use GROOT_Global_m, only: nT
    use GROOT_Parameters_m, only: clen
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    type(orion_data),        intent(in)    :: IOfield_wall
    integer :: ib, ifa, cc, Nxc, Nyc, Nzc
    integer :: n1, n2, i1, i2, i, j, k, facetype
    character(len=clen) :: facename
    character(1)        :: bname, fname
    character(100)      :: subargs(2)

    if (allocated(viscous_flag)) deallocate(viscous_flag)
    allocate(viscous_flag(1:grid%nb, 1:6))
    viscous_flag = .false.

    do ib = 1, grid%nb
      Nxc = grid%blk(ib)%dim(1)
      Nyc = grid%blk(ib)%dim(2)
      Nzc = grid%blk(ib)%dim(3)

      !! Override wall temperatures from IOfield_wall
      do ifa = 1, 6
        write(bname, '(I1)') ib
        write(fname, '(I1)') ifa
        facename = 'B'//bname//'F'//fname

        do cc = 1, size(IOfield_wall%block)
          if (trim(facename) /= IOfield_wall%block(cc)%name) cycle
          viscous_flag(ib, ifa) = .true.
          write(*, '(A,A)') '   Found viscous face ', facename

          call parse(facename, 'F', subargs)
          read(subargs(2), '(I1)') facetype
          select case (facetype)
            case (1,2);  n1 = Nyc; n2 = Nzc
            case (3,4);  n1 = Nxc; n2 = Nzc
            case (5,6);  n1 = Nxc; n2 = Nyc
          end select

          do i2 = 1, n2; do i1 = 1, n1
            select case (ifa)
              case (1,2);  i=1;  j=i1; k=i2
              case (3,4);  i=i1; j=1;  k=i2
              case (5,6);  i=i1; j=i2; k=1
            end select
            if (grid%blk(ib)%face(ifa)%bc(i1,i2) > 4) then
              grid%blk(ib)%face(ifa)%T(i1,i2) = IOfield_wall%block(cc)%vars(nT, i, j, k)
            else if (grid%blk(ib)%face(ifa)%bc(i1,i2) == 1) then
              grid%blk(ib)%face(ifa)%T(i1,i2) = 0.0_R8
            end if
          end do; end do
        end do
      end do
    end do
  end subroutine Setup_Wall


  !! ----------------------------------------------------------
  !! Internal helper: fill face Q, Qtot, eps, a, T defaults.
  subroutine Setup_Face(ib, ifa, nc, grid)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use GROOT_Global_m, only: eps_wall, Ngg
    implicit none
    integer,                 intent(in)    :: ib, ifa
    integer,                 intent(in)    :: nc(2)
    type(GROOT_domain_type), intent(inout) :: grid
    integer :: n1, n2

    n1 = nc(1); n2 = nc(2)
    grid%blk(ib)%face(ifa)%n = nc
    grid%blk(ib)%face(ifa)%Q   (1:n1, 1:n2, 0:Ngg) = 0.0_R8
    grid%blk(ib)%face(ifa)%Qtot(1:n1, 1:n2)         = 0.0_R8
    grid%blk(ib)%face(ifa)%eps (1:n1, 1:n2)         = eps_wall
    grid%blk(ib)%face(ifa)%a   (1:n1, 1:n2, 0:Ngg) = 1.0_R8
    grid%blk(ib)%face(ifa)%T   (1:n1, 1:n2)         = 0.0_R8
  end subroutine Setup_Face


  !! ----------------------------------------------------------
  !! Update face%a with Planck fractions at wall temperature for SNBW.
  !!
  !! Setup_Face initialises face%a = 1 for all bands.  For WSGG (Ngg=4)
  !! the error is tolerable, but for SNBW (Ngg=450) a=1 for all bands
  !! would overestimate wall emission by a factor of 450.  This routine
  !! computes the correct Planck fractions a(ig,T_wall) for each wall face.
  !! Non-wall faces (T=0 from initialisation) are skipped.
  subroutine Setup_SNBWall(grid)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use GROOT_Lib_SNBW,         only: snb_planck_weight, N_SNB_BANDS
    use GROOT_Global_m,         only: Ngg
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid

    integer  :: ib, ifa, i1, i2, ig
    real(R8) :: T

    do ib = 1, grid%nb
      do ifa = 1, 6
        do i2 = 1, grid%blk(ib)%face(ifa)%n(2)
          do i1 = 1, grid%blk(ib)%face(ifa)%n(1)
            T = grid%blk(ib)%face(ifa)%T(i1, i2)
            if (T < 1.0_R8) cycle   ! non-wall or uninitialised face
            grid%blk(ib)%face(ifa)%a(i1, i2, 0) = 0.0_R8
            do ig = 1, Ngg
              grid%blk(ib)%face(ifa)%a(i1, i2, ig) = snb_planck_weight(T, ig)
            end do
          end do
        end do
      end do
    end do
  end subroutine Setup_SNBWall

  !! ----------------------------------------------------------
  !! Update face%a with WSGG polynomial weights at wall temperature.
  !!
  !! Analogous to Setup_SNBWall for SNB/SNBW. Without this, Setup_Face leaves
  !! face%a = 1.0 for all bands, so each band emits the full blackbody power
  !! instead of the correct fraction a_wsgg(ig,T_wall). The total wall emission
  !! is overestimated by a factor of (Ngg+1) = 5.
  !!
  !! For 'wsgg-H2O':    a(ig) depends only on T (not xH2O or p).
  !! For 'wsgg-H2OCO2': a(ig) also depends on MR=xH2O/xCO2; MR=2 is used as
  !!                    representative (typical for H2O-CO2 combustion products).
  subroutine Setup_WSGGWall(grid)
    use GROOT_Advanced_Types_m,  only: GROOT_domain_type
    use GROOT_Lib_Radproperties, only: wsggH2O, wsggH2OCO2_MR
    use GROOT_Global_m,          only: Ngg, model
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid

    integer  :: ib, ifa, i1, i2
    real(R8) :: T_wall, kp_dum(0:4), a_w(0:4), eps_dum
    ! MR=2 reference: xH2O = 2/3, xCO2 = 1/3
    real(R8), parameter :: XH2O_REF = 2.0_R8 / 3.0_R8
    real(R8), parameter :: XCO2_REF = 1.0_R8 / 3.0_R8

    do ib = 1, grid%nb
      do ifa = 1, 6
        do i2 = 1, grid%blk(ib)%face(ifa)%n(2)
          do i1 = 1, grid%blk(ib)%face(ifa)%n(1)
            T_wall = grid%blk(ib)%face(ifa)%T(i1, i2)
            if (T_wall < 1.0_R8) cycle   ! non-wall or uninitialised face

            if (trim(model) == 'wsgg-H2O') then
              ! a(ig) depends only on T; xH2O and p are dummies here
              call wsggH2O(1.0_R8, 1.0e5_R8, T_wall, 1.0_R8, kp_dum, a_w, eps_dum)
            elseif (trim(model) == 'wsgg-H2OCO2') then
              ! a(ig) depends on T and MR=xH2O/xCO2; MR=2 representative
              call wsggH2OCO2_MR(XH2O_REF, XCO2_REF, 1.0e5_R8, T_wall, 1.0_R8, &
                                  kp_dum, a_w, eps_dum)
            else
              cycle
            end if
            grid%blk(ib)%face(ifa)%a(i1, i2, 0:Ngg) = a_w(0:Ngg)
          end do
        end do
      end do
    end do
  end subroutine Setup_WSGGWall


  !! ----------------------------------------------------------
  !! Convert partial densities to molar fractions in-place
  subroutine Compute_Molfractions(grid)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use GROOT_Global_m, only: nsc, Ru_, wm_tab
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    integer  :: b, i, j, k, is
    real(R8) :: rhotot, wm

    do b = 1, grid%nb
      associate(Nx => grid%blk(b)%dim(1), &
                Ny => grid%blk(b)%dim(2), &
                Nz => grid%blk(b)%dim(3))
        do k = 1, Nz; do j = 1, Ny; do i = 1, Nx
          rhotot = sum(grid%blk(b)%gas_phase%rho(1:nsc, i, j, k))
          wm     = Ru_ / grid%blk(b)%gas_phase%R(i, j, k)
          do is = 1, nsc
            grid%blk(b)%gas_phase%rho(is,i,j,k) = &
              grid%blk(b)%gas_phase%rho(is,i,j,k) / rhotot * wm / wm_tab(is)
          end do
        end do; end do; end do
      end associate
    end do
  end subroutine Compute_Molfractions

end module GROOT_Mod_Phase
