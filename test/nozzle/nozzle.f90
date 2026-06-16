!! ============================================================
!! Nozzle wall heat flux test.
!!
!! Reads a 1-D axial profile from nozzle.dat (145 stations with
!! X, R_wall, T, p[bar], xH2O, xCO2 from CEA), builds a 2-D
!! axisymmetric structured mesh (Ni_ax=50, Nj_rad=20, Nk=1 wedge),
!! and computes the lateral wall heat flux with four spectral models:
!!   gray         — gray-gas Planck-mean absorption (H2O + CO2)
!!   wsgg-H2OCO2  — WSGG by Fabiani et al. (JPP 2025), 4 gray gases
!!   snbw         — Statistical Narrow Band, weak (Beer-Lambert)
!!   snb          — Statistical Narrow Band, full Malkmus
!!
!! Gas properties vary in the axial direction only.
!! No CO in the mixture (xCO = 0 throughout).
!!
!! Wall conditions: cold black walls everywhere (T_wall=0 K, eps=1)
!!   → Qtot on face 4 (lateral) is the incident radiative flux.
!!
!! Output: OUTPUT/heatflux_nozzle_<model>.dat
!!   Columns: x[m]  Qtot[W/m2]
!! ============================================================
program nozzle_test
  use iso_fortran_env, only: I4 => int32, R8 => real64
#if defined (_OPENMP)
  use omp_lib
#endif
  use GROOT_Config_Types_m,    only: obj_dtm, obj_rad_model, obj_io, obj_sim_param
  use GROOT_Assign_Setup,      only: Assign_Setup
  use GROOT_Advanced_Types_m,  only: GROOT_simulation_type
  use GROOT_Global_m
  use GROOT_Mod_Allocate_Data, only: Setup_Data_Structure
  use GROOT_Mod_Metrics,       only: Setup_Metrics
  use GROOT_Wrap_Solve,        only: GROOT_solve
  use GROOT_Lib_Radproperties, only: compute_gasproperties, map_species
  use GROOT_Mod_Phase,         only: Setup_WSGGWall, Setup_SNBWall
  use Lib_ORION_Data
  implicit none

  !! ── Mesh parameters ────────────────────────────────────────
  integer,  parameter :: Ni_ax  = 100    !< axial cells
  integer,  parameter :: Nj_rad = 50    !< radial cells
  real(R8), parameter :: dth    = 0.01_R8  !< full wedge angle [rad]
  real(R8), parameter :: pi_val = acos(-1.0_R8)

  !! ── Nozzle data ─────────────────────────────────────────────
  integer,          parameter :: N_DAT   = 145

  !! ── Profile arrays ──────────────────────────────────────────
  real(R8) :: x_dat  (N_DAT), R_dat  (N_DAT)
  real(R8) :: T_dat  (N_DAT), p_bar_dat(N_DAT)
  real(R8) :: xH2O_dat(N_DAT), xCO2_dat(N_DAT)

  !! ── Interpolated mesh arrays ─────────────────────────────────
  real(R8) :: x_node(0:Ni_ax), R_node(0:Ni_ax)
  real(R8) :: x_cell (Ni_ax), R_cell (Ni_ax)
  real(R8) :: T_cell (Ni_ax), p_Pa_cell(Ni_ax)
  real(R8) :: xH2O_cell(Ni_ax), xCO2_cell(Ni_ax)

  type(GROOT_simulation_type) :: simulation
  type(orion_data)             :: IOgas
  real(R8) :: t1, t2
  integer  :: im

  integer,  parameter  :: N_MODELS = 4
  character(len=16) :: model_list(N_MODELS)
  model_list(1) = 'gray'
  model_list(2) = 'wsgg-H2OCO2'
  model_list(3) = 'snbw'
  model_list(4) = 'snb'

  !! ── OpenMP ──────────────────────────────────────────────────
#if defined (_OPENMP)
  !$omp parallel
  obj_sim_param%nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  write(*,'(A,I4,A)') ' OpenMP parallel: ', obj_sim_param%nthreads, ' threads'
#else
  obj_sim_param%nthreads = 1
  write(*,'(A)') ' Serial execution'
#endif

  !! ── Species globals (H2O + CO2, no CO) ─────────────────────
  !! These are shared across all model runs; only nscdat / ind_species
  !! change between gray and the other models (see run_model).
  nsc = 2; nsoot = 0; nrans = 0; nscrad = 0
  allocate(cfd_species(1:nsc), wm_tab(1:nsc))
  cfd_species(1) = 'H2O';  wm_tab(1) = 18.015_R8
  cfd_species(2) = 'CO2';  wm_tab(2) = 44.010_R8

  !! ── Common DTM options ──────────────────────────────────────
  obj_dtm%Nr              = 400
  obj_dtm%tol             = 1.0e-3_R8
  obj_dtm%itermax         = 1          !< black walls → one-shot solve
  obj_dtm%source          = .false.    !< wall flux only, skip volumetric source
  obj_dtm%twoDax          = .true.
  obj_dtm%wedge_angle_deg = dth * 180.0_R8 / pi_val
  obj_dtm%optically_thin  = .false.
  obj_io%sol_format       = 'tecplot ascii'
  obj_io%gas_format       = 'tecplot ascii'
  obj_rad_model%eps_wall  = 1.0_R8    !< eps, overridden per-face in setup_bc

  !! ── Load profile ────────────────────────────────────────────
  write(*,'(/,A)') ' Reading nozzle.dat ...'
  call read_profile(x_dat, R_dat, T_dat, p_bar_dat, xH2O_dat, xCO2_dat)

  !! ── Interpolate to 50 axial cells ───────────────────────────
  call interpolate_profile(x_dat, R_dat, T_dat, p_bar_dat, xH2O_dat, xCO2_dat, &
                            x_node, R_node, x_cell, R_cell,                       &
                            T_cell, p_Pa_cell, xH2O_cell, xCO2_cell)

  !! ── Build mesh (geometry only, identical for all models) ─────
  call build_iogas(IOgas, x_node, R_node)

  !! ── Output directory ─────────────────────────────────────────
  call execute_command_line('mkdir -p OUTPUT', wait=.true.)

  write(*,'(/,A)') ' ================================================='
  write(*,'(  A)') ' Nozzle wall heat flux — 4 spectral models'
  write(*,'(  A,I0,A,I0,A)') ' Mesh: ', Ni_ax, ' x ', Nj_rad, ' x 1 (2-D axi)'
  write(*,'(  A,F8.4,A,F8.4,A)') ' X range: [', x_dat(1), ',', x_dat(N_DAT), '] m'
  write(*,'(  A)') ' ================================================='

  !! ── Run each spectral model ──────────────────────────────────
  do im = 1, N_MODELS
    write(*,'(/,A,A)') ' ─── Model: ', trim(model_list(im))
    call run_model(trim(adjustl(model_list(im))), simulation, IOgas, &
                   x_cell, T_cell, p_Pa_cell, xH2O_cell, xCO2_cell, t1, t2)
    write(*,'(A,F8.4,A)') '   CPU time: ', &
      (t2-t1)/real(obj_sim_param%nthreads,R8)/60.0_R8, ' min'
    call write_output(simulation%domain, trim(adjustl(model_list(im))), x_cell)
  end do

  write(*,'(/,A)') ' Done. Results in OUTPUT/heatflux_nozzle_<model>.dat'

contains

  !! ============================================================
  !! Read nozzle profile from dat file.
  !! Skip VARIABLES header line; read N_DAT rows, columns 1–6.
  !! ============================================================
  subroutine read_profile(x, R, T, p_bar, xH2O, xCO2)
    real(R8), intent(out) :: x(N_DAT), R(N_DAT), T(N_DAT)
    real(R8), intent(out) :: p_bar(N_DAT), xH2O(N_DAT), xCO2(N_DAT)
    real(R8) :: dum(7)
    integer  :: u, i

    character(len=64) :: dat_path
    logical :: found
    dat_path = 'nozzle.dat'
    inquire(file=trim(dat_path), exist=found)
    if (.not. found) then
      dat_path = '../nozzle.dat'
      inquire(file=trim(dat_path), exist=found)
    end if
    if (.not. found) then
      write(*,'(A)') ' ERROR: nozzle.dat not found in ./ or ../'
      stop 1
    end if
    open(newunit=u, file=trim(dat_path), status='old', action='read')
    read(u, *)        ! skip VARIABLES header
    do i = 1, N_DAT
      read(u,*) x(i), R(i), T(i), p_bar(i), xH2O(i), xCO2(i), &
                dum(1), dum(2), dum(3), dum(4), dum(5), dum(6), dum(7)
    end do
    close(u)

    write(*,'(A,F8.4,A,F8.4,A)') '   x    range: [', x(1),      ',', x(N_DAT),      '] m'
    write(*,'(A,F8.4,A,F8.4,A)') '   R    range: [', minval(R),  ',', maxval(R),      '] m'
    write(*,'(A,F6.1,A,F6.1,A)') '   T    range: [', minval(T),  ',', maxval(T),      '] K'
    write(*,'(A,F8.4,A,F8.4,A)') '   p    range: [', minval(p_bar), ',', maxval(p_bar), '] bar'
    write(*,'(A,F6.4,A,F6.4,A)') '   xH2O range: [', minval(xH2O),',', maxval(xH2O), ']'
    write(*,'(A,F6.4,A,F6.4,A)') '   xCO2 range: [', minval(xCO2),',', maxval(xCO2), ']'
  end subroutine read_profile


  !! ============================================================
  !! Piecewise-linear interpolation of scalar y_d(x_d) at xq.
  !! ============================================================
  function lininterp(x_d, y_d, xq) result(yq)
    real(R8), intent(in) :: x_d(N_DAT), y_d(N_DAT), xq
    real(R8) :: yq, t
    integer  :: k
    if (xq <= x_d(1))     then; yq = y_d(1);     return; end if
    if (xq >= x_d(N_DAT)) then; yq = y_d(N_DAT); return; end if
    do k = 1, N_DAT - 1
      if (xq >= x_d(k) .and. xq <= x_d(k+1)) then
        t  = (xq - x_d(k)) / (x_d(k+1) - x_d(k))
        yq = y_d(k) + t * (y_d(k+1) - y_d(k))
        return
      end if
    end do
    yq = y_d(N_DAT)
  end function lininterp


  !! ============================================================
  !! Interpolate profile arrays to evenly-spaced nodes and cell centres.
  !! ============================================================
  subroutine interpolate_profile(x_d, R_d, T_d, p_d, xH2O_d, xCO2_d, &
                                  x_n, R_n, x_c, R_c, T_c, p_Pa_c,    &
                                  xH2O_c, xCO2_c)
    real(R8), intent(in)  :: x_d(N_DAT), R_d(N_DAT), T_d(N_DAT)
    real(R8), intent(in)  :: p_d(N_DAT), xH2O_d(N_DAT), xCO2_d(N_DAT)
    real(R8), intent(out) :: x_n(0:Ni_ax), R_n(0:Ni_ax)
    real(R8), intent(out) :: x_c(Ni_ax), R_c(Ni_ax), T_c(Ni_ax)
    real(R8), intent(out) :: p_Pa_c(Ni_ax), xH2O_c(Ni_ax), xCO2_c(Ni_ax)
    real(R8) :: dx
    integer  :: i

    dx = (x_d(N_DAT) - x_d(1)) / real(Ni_ax, R8)

    do i = 0, Ni_ax
      x_n(i) = x_d(1) + i * dx
      R_n(i) = lininterp(x_d, R_d, x_n(i))
    end do

    do i = 1, Ni_ax
      x_c(i)    = x_d(1) + (i - 0.5_R8) * dx
      R_c(i)    = lininterp(x_d, R_d,    x_c(i))
      T_c(i)    = lininterp(x_d, T_d,    x_c(i))
      p_Pa_c(i) = lininterp(x_d, p_d,    x_c(i)) * 1.0e5_R8  ! bar → Pa
      xH2O_c(i) = lininterp(x_d, xH2O_d, x_c(i))
      xCO2_c(i) = lininterp(x_d, xCO2_d, x_c(i))
    end do
  end subroutine interpolate_profile


  !! ============================================================
  !! Build IOgas: nozzle wedge mesh (geometry only, no field vars).
  !! i-direction : axial   (x),    nodes at x_n(0:Ni_ax)
  !! j-direction : radial  (r),    0 at axis → R_n(i) at wall
  !! k-direction : azimuth (wedge), Nk=1
  !!
  !! Node coords:
  !!   mesh(1) = x_n(i)
  !!   mesh(2) = r * cos(dth/2)
  !!   mesh(3) = ±r * sin(dth/2)   (k=0 → -, k=1 → +)
  !! ============================================================
  subroutine build_iogas(IOg, x_n, R_n)
    type(orion_data), intent(out) :: IOg
    real(R8),         intent(in)  :: x_n(0:Ni_ax), R_n(0:Ni_ax)
    integer  :: i, j, k
    real(R8) :: rr, cdth, sdth
    real(R8), parameter :: r_min = 1.0e-6_R8  ! avoid degenerate axis cells

    cdth = cos(dth * 0.5_R8)
    sdth = sin(dth * 0.5_R8)

    allocate(IOg%block(1))
    IOg%block(1)%Ni   = Ni_ax
    IOg%block(1)%Nj   = Nj_rad
    IOg%block(1)%Nk   = 1
    IOg%block(1)%name = 'Block1'

    allocate(IOg%block(1)%mesh(1:3, 0:Ni_ax, 0:Nj_rad, 0:1))

    do k = 0, 1; do j = 0, Nj_rad; do i = 0, Ni_ax
      rr = max(real(j, R8) * R_n(i) / real(Nj_rad, R8), r_min)
      IOg%block(1)%mesh(1, i, j, k) = x_n(i)
      IOg%block(1)%mesh(2, i, j, k) = rr * cdth
      IOg%block(1)%mesh(3, i, j, k) = rr * sdth * real(2*k - 1, R8)
    end do; end do; end do
  end subroutine build_iogas


  !! ============================================================
  !! Run one spectral model: configure → mesh → BCs → rad props → solve.
  !! ============================================================
  subroutine run_model(mname, sim, IOg, x_c, T_c, p_Pa_c, xH2O_c, xCO2_c, &
                        t_start, t_end)
    character(len=*),            intent(in)    :: mname
    type(GROOT_simulation_type), intent(inout) :: sim
    type(orion_data),            intent(inout) :: IOg
    real(R8),                    intent(in)    :: x_c(Ni_ax)
    real(R8),                    intent(in)    :: T_c(Ni_ax), p_Pa_c(Ni_ax)
    real(R8),                    intent(in)    :: xH2O_c(Ni_ax), xCO2_c(Ni_ax)
    real(R8),                    intent(out)   :: t_start, t_end

    !! ── 1. Configure model and sync globals ──
    obj_rad_model%model  = mname
    obj_rad_model%k_user = 1.0_R8   ! irrelevant for non-const models
    call Assign_Setup()              ! sets Ngg, model, Nr, twoDax, …

    !! ── 2. Species: per-model handling ──
    if (allocated(rad_species)) deallocate(rad_species)
    if (allocated(ind_species)) deallocate(ind_species)
    select case (trim(mname))
      case ('gray')
        !! graygas(T,p,1) = H2O,  graygas(T,p,2) = CO2
        nscdat = 2
        allocate(rad_species(1:nscdat), ind_species(1:nscdat))
        rad_species(1) = 'H2O';  ind_species(1) = 1
        rad_species(2) = 'CO2';  ind_species(2) = 2
      case default
        !! WSGG/SNB: map_species sets module-level iH2O, iCO2, iCO
        !! iCO = findindex('CO') = 0 because CO is not in cfd_species → xCO = 0
        nscdat = 0
        allocate(rad_species(0), ind_species(0))
        call map_species()
    end select

    !! ── 3. Allocate domain (deallocates previous run automatically) ──
    call Setup_Data_Structure(sim%domain, IOg)
    call Setup_Metrics(sim%domain, IOg)

    !! ── 4. Boundary conditions ──
    call setup_bc_nozzle(sim%domain)

    !! ── 5. Radiative properties ──
    call setup_rad_nozzle(sim%domain, T_c, p_Pa_c, xH2O_c, xCO2_c)

    !! ── 6. Wall spectral weights (fix face%a from default 1.0) ──
    select case (trim(mname))
      case ('wsgg-H2O', 'wsgg-H2OCO2')
        call Setup_WSGGWall(sim%domain)
      case ('snbw', 'snb')
        call Setup_SNBWall(sim%domain)
    end select

    !! ── 7. Solve ──
    call cpu_time(t_start)
    call GROOT_solve(sim)
    call cpu_time(t_end)
  end subroutine run_model


  !! ============================================================
  !! Nozzle boundary conditions (same for all models):
  !!   Face 1 (x=x_min) : inlet endcap  — cold black wall
  !!   Face 2 (x=x_max) : outlet endcap — cold black wall
  !!   Face 3 (r=0)      : axis          — specular (bc=300)
  !!   Face 4 (r=R_wall) : lateral wall  — cold black wall
  !!   Face 5,6 (azimuth): periodic      — bc=200
  !!
  !! Cold walls (T=0) → Qtot on face 4 = incident radiative flux.
  !! ============================================================
  subroutine setup_bc_nozzle(grid)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    type(GROOT_domain_type), intent(inout) :: grid
    integer :: b

    do b = 1, grid%nb
      !! Inlet
      grid%blk(b)%face(1)%bc (:,:) = 301
      grid%blk(b)%face(1)%T  (:,:) = 0.0_R8
      grid%blk(b)%face(1)%eps(:,:) = 1.0_R8

      !! Outlet
      grid%blk(b)%face(2)%bc (:,:) = 301
      grid%blk(b)%face(2)%T  (:,:) = 0.0_R8
      grid%blk(b)%face(2)%eps(:,:) = 1.0_R8

      !! Axis (specular reflection)
      grid%blk(b)%face(3)%bc (:,:) = 300
      grid%blk(b)%face(3)%T  (:,:) = 0.0_R8
      grid%blk(b)%face(3)%eps(:,:) = 0.0_R8

      !! Lateral wall — cold black (incident flux = Qtot)
      grid%blk(b)%face(4)%bc (:,:) = 301
      grid%blk(b)%face(4)%T  (:,:) = 0.0_R8
      grid%blk(b)%face(4)%eps(:,:) = 1.0_R8

      !! Azimuthal (wedge faces — periodic / axisymmetric)
      grid%blk(b)%face(5)%bc (:,:) = 200
      grid%blk(b)%face(5)%T  (:,:) = 0.0_R8
      grid%blk(b)%face(5)%eps(:,:) = 0.0_R8

      grid%blk(b)%face(6)%bc (:,:) = 200
      grid%blk(b)%face(6)%T  (:,:) = 0.0_R8
      grid%blk(b)%face(6)%eps(:,:) = 0.0_R8
    end do
  end subroutine setup_bc_nozzle


  !! ============================================================
  !! Fill radiative property arrays per cell.
  !! Properties vary in the axial direction (i) only; all radial
  !! stations (j) at the same axial position share identical values.
  !! ============================================================
  subroutine setup_rad_nozzle(grid, T_c, p_Pa_c, xH2O_c, xCO2_c)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    type(GROOT_domain_type), intent(inout) :: grid
    real(R8), intent(in) :: T_c(Ni_ax), p_Pa_c(Ni_ax)
    real(R8), intent(in) :: xH2O_c(Ni_ax), xCO2_c(Ni_ax)

    real(R8) :: xvec(nsc), ka_loc(0:Ngg), a_loc(0:Ngg), Ib_loc
    real(R8) :: kb_loc(0:Ngg)
    logical  :: snb_full
    integer  :: b, i, j, Nx, Ny, Nz

    snb_full = (trim(model) == 'snb')

    do b = 1, grid%nb
      Nx = grid%blk(b)%dim(1)
      Ny = grid%blk(b)%dim(2)
      Nz = grid%blk(b)%dim(3)

      allocate(grid%blk(b)%ka    (1:Nx, 1:Ny, 1:Nz, 0:Ngg))
      allocate(grid%blk(b)%a     (1:Nx, 1:Ny, 1:Nz, 0:Ngg))
      allocate(grid%blk(b)%Ib   (1:Nx, 1:Ny, 1:Nz))
      allocate(grid%blk(b)%source(1:Nx, 1:Ny, 1:Nz))
      allocate(grid%blk(b)%G     (1:Nx, 1:Ny, 1:Nz))
      if (snb_full) allocate(grid%blk(b)%kb(1:Nx, 1:Ny, 1:Nz, 0:Ngg))

      grid%blk(b)%source = 0.0_R8
      grid%blk(b)%G      = 0.0_R8

      do i = 1, Nx
        xvec(1) = xH2O_c(i)   ! molar fraction H2O
        xvec(2) = xCO2_c(i)   ! molar fraction CO2

        if (snb_full) then
          call compute_gasproperties(xvec, p_Pa_c(i), T_c(i), 0.0_R8, &
                                      ka_loc, a_loc, Ib_loc, kb=kb_loc)
        else
          call compute_gasproperties(xvec, p_Pa_c(i), T_c(i), 0.0_R8, &
                                      ka_loc, a_loc, Ib_loc)
        end if

        do j = 1, Ny
          grid%blk(b)%ka (i, j, 1, 0:Ngg) = ka_loc
          grid%blk(b)%a  (i, j, 1, 0:Ngg) = a_loc
          grid%blk(b)%Ib (i, j, 1)         = Ib_loc
          if (snb_full) grid%blk(b)%kb(i, j, 1, 0:Ngg) = kb_loc
        end do
      end do
    end do
  end subroutine setup_rad_nozzle


  !! ============================================================
  !! Write lateral wall heat flux (face 4) vs axial position.
  !! ============================================================
  subroutine write_output(grid, mname, x_c)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    type(GROOT_domain_type), intent(in) :: grid
    character(len=*),        intent(in) :: mname
    real(R8),                intent(in) :: x_c(Ni_ax)
    integer :: u, i

    open(newunit=u, file='OUTPUT/heatflux_nozzle_'//trim(mname)//'.dat', &
         action='write')
    write(u,'(A)') '# Nozzle lateral wall heat flux (face 4, r=R_wall)'
    write(u,'(A)') '# Model: '//trim(mname)
    write(u,'(A)') '#   x[m]              Qtot[W/m2]'
    do i = 1, Ni_ax
      write(u,'(2ES20.8)') x_c(i), grid%blk(1)%face(4)%Qtot(i, 1)
    end do
    close(u)

    write(*,'(A)') '   → OUTPUT/heatflux_nozzle_'//trim(mname)//'.dat'
  end subroutine write_output

end program nozzle_test
