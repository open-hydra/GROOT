!! ============================================================
!! Axisymmetric furnace test with non-uniform temperature field.
!!
!! Geometry : radius R=0.45 m, axial length L=5 m.
!!            Thin wedge: δθ=0.01 rad, Nk=1 (axisymmetric).
!!            Ni_ax x Nj_rad x 1 cells.
!!
!! Medium   : non-uniform T(x,r) from Tfun (polynomial profile
!!            fitted to a turbulent flame, Tref≈1563 K at peak).
!!            Absorption coeff k = k_user (const model).
!!
!! BCs      :
!!   Face 1 (x=0)    → black endcap  (bc=301, Tw1=820.11 K)
!!   Face 2 (x=L)    → black endcap  (bc=301, Tw2=615.19 K)
!!   Face 3 (r=0)    → axis          (bc=300, specular)
!!   Face 4 (r=R)    → black lateral (bc=301, Tw_lat=349.49 K)
!!   Face 5,6 (azim) → periodic      (bc=200)
!!
!! Outputs (vs axial position x):
!!   OUTPUT/source_j25.dat   source at j≈25%R (r≈0.101 m)
!!   OUTPUT/source_j75.dat   source at j≈75%R (r≈0.326 m)
!!   OUTPUT/heatflux_lat.dat lateral wall heat flux on face 4
!!
!! Reference : data/ref_source25.dat, ref_source75.dat, ref_sourceqw.dat
!! ============================================================
program furnace_test
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
  use Lib_ORION_Data

  implicit none

  !! Geometry and mesh parameters
  integer,  parameter :: Ni_ax  = 200          !< axial  cells  (matches old 200x100 run)
  integer,  parameter :: Nj_rad = 100          !< radial cells
  real(R8), parameter :: R_cyl  = 0.45_R8      !< cylinder radius [m]
  real(R8), parameter :: L_cyl  = 5.0_R8       !< cylinder length [m]
  real(R8), parameter :: dth    = 0.01_R8      !< wedge angle     [rad]

  !! Wall temperatures [K]
  real(R8), parameter :: Tw_face1 = 820.1117_R8   !< x=0 endcap
  real(R8), parameter :: Tw_face2 = 615.1889_R8   !< x=L endcap
  real(R8), parameter :: Tw_lat   = 349.4948_R8   !< lateral wall

  !! Radial output indices (≈25% and ≈75% of radius)
  integer,  parameter :: j_25 = 25   !< r_c = 24.5/100*R = 0.110 m = 24.5% R
  integer,  parameter :: j_75 = 75   !< r_c = 74.5/100*R = 0.335 m = 74.5% R

  type(GROOT_simulation_type) :: simulation
  type(orion_data)             :: IOgas
  real(R8) :: t1, t2

  !! --------------------------------------------------------
  !! Thread detection
#if defined (_OPENMP)
  !$omp parallel
  obj_sim_param%nthreads = OMP_GET_NUM_THREADS()
  !$omp end parallel
  write(*,'(A)')    ' Parallel execution (OpenMP)'
  write(*,'(A,I4)') '   Threads: ', obj_sim_param%nthreads
#else
  obj_sim_param%nthreads = 1
  write(*,'(A)') ' Serial execution'
#endif

  !! --------------------------------------------------------
  !! Configure (2-D axisymmetric DTM, const model)
  obj_dtm%Nr             = 256
  obj_dtm%tol            = 1.0e-3_R8
  obj_dtm%itermax        = 1
  obj_dtm%source         = .true.
  obj_dtm%twoDax         = .true.
  obj_dtm%wedge_angle_deg= dth * 180.0_R8 / (4.0_R8 * datan(1.0_R8))
  obj_dtm%optically_thin = .false.

  obj_rad_model%model    = 'const'
  obj_rad_model%eps_wall = 1.0_R8
  obj_rad_model%k_user   = 0.3471_R8

  obj_io%sol_format = 'tecplot ascii'
  obj_io%gas_format = 'tecplot ascii'


  !! Set species manually
  nsc  = 1
  nsoot= 0
  nrans = 0
  allocate(cfd_species(1:nsc), wm_tab(1:nsc))
  cfd_species(1) = 'N2'
  wm_tab(1)      = 28.0_R8
  nscrad = 0
  nscdat = 0
  allocate(rad_species(0), ind_species(0))

  call Assign_Setup()

  write(*,'(A)')
  write(*,'(A)') ' =============================================='
  write(*,'(A)') ' Axisymmetric furnace test (non-uniform T)'
  write(*,'(A)') '   Domain  : R=0.45 m, L=5 m'
  write(*,'(A,I0,A,I0,A)') '   Mesh    : ', Ni_ax, ' x ', Nj_rad, ' x 1  (same as ref)'
  write(*,'(A,F10.6,A)') '   Wedge   : ', dth, ' rad'
  write(*,'(A)') '   T(x,r)  : polynomial furnace profile (~1563 K peak)'
  write(*,'(A,F8.4)') '   k_user  = ', k_user
  write(*,'(A,I0)') '   Nr      = ', Nr
  write(*,'(A)') ' =============================================='
  write(*,'(A)')

  !! --------------------------------------------------------
  !! Build IOgas in memory: wedge mesh + non-uniform T field
  call Create_IOgas(IOgas, Ni_ax, Nj_rad, R_cyl, L_cyl)

  !! --------------------------------------------------------
  !! Allocate grid and compute mesh metrics
  call Setup_Data_Structure(simulation%domain, IOgas)
  call Setup_Metrics(simulation%domain, IOgas)

  !! --------------------------------------------------------
  !! Set BCs (different Tw on each face)
  call Setup_BC(simulation%domain)

  !! --------------------------------------------------------
  !! Set radiative properties (non-uniform Ib from Tfun)
  call Setup_Rad_NonUniform(simulation%domain, Ni_ax, Nj_rad, R_cyl, L_cyl)

  !! --------------------------------------------------------
  !! Solve
  call cpu_time(t1)
  call GROOT_solve(simulation)
  call cpu_time(t2)

  write(*,'(A)')
  write(*,'(A,F8.4,A)') '   CPU time: ', &
    (t2-t1)/real(obj_sim_param%nthreads,R8)/60.0_R8, ' min'

  !! --------------------------------------------------------
  call Write_Output(simulation%domain, Ni_ax, L_cyl)

  write(*,'(A)') ' Done. See OUTPUT/source_j25.dat, source_j75.dat, heatflux_lat.dat'
  write(*,'(A)')

contains

  !! ----------------------------------------------------------
  !! Temperature profile from the furnace benchmark
  !! (polynomial fit, Tref=1563 K, R=0.45 m, L=5 m).
  !! Arguments: x_ax = axial coordinate, r_rad = radial coordinate.
  pure function Tfun(x_ax, r_rad) result(T)
    implicit none
    real(R8), intent(in) :: x_ax, r_rad
    real(R8) :: T
    real(R8) :: Rf, Lf, zpmax, Lz, L0, ztilde, zp
    real(R8) :: Tref, Ti, Te, Tr, de, f, a
    Rf    = 0.45_R8
    Lf    = 5.0_R8
    zpmax = 0.64_R8
    Lz    = 5.5555_R8
    L0    = 0.45_R8
    Tref  = 1563.0_R8
    Ti    = 0.1906_R8
    Te    = 0.6769_R8
    Tr    = 0.6865_R8
    de    = -0.2924_R8

    ztilde = (x_ax - Lf * 0.5_R8) / L0
    zp     = ztilde / Lz
    f      = 1.0_R8 - 3.0_R8*(r_rad/Rf)**2 + 2.0_R8*(r_rad/Rf)**3

    if (zp <= -zpmax) then
      a = 1.0_R8 + (1.0_R8 - Ti) * ((zp + zpmax) / (1.0_R8 - zpmax))**3
    else
      a = 1.0_R8 &
        - (de*(1.0_R8+zpmax) + 3.0_R8*(1.0_R8-Te)) * ((zp+zpmax)/(1.0_R8+zpmax))**2 &
        + (de*(1.0_R8+zpmax) + 2.0_R8*(1.0_R8-Te)) * ((zp+zpmax)/(1.0_R8+zpmax))**3
    end if
    T = ((a - Tr) * f + Tr) * Tref
  end function Tfun


  !! ----------------------------------------------------------
  !! Build axisymmetric wedge cylinder mesh with non-uniform T.
  !!   i-direction : axial   (x), 0 to L_cyl
  !!   j-direction : radial  (y projected), 0 to R_cyl
  !!   k-direction : azimuth (z projected), Nk=1
  !!   mesh(2) = r * cos(dth/2),  mesh(3) = ±r * sin(dth/2)
  subroutine Create_IOgas(IOgas, Ni, Nj, R0, L0)
    implicit none
    type(orion_data), intent(out) :: IOgas
    integer,  intent(in) :: Ni, Nj
    real(R8), intent(in) :: R0, L0
    integer  :: i, j, k
    real(R8) :: p0, R_air, rho, rr, cdth, sdth, x_c, r_c
    real(R8), parameter :: r_min = 1.0e-6_R8
    real(R8), parameter :: T_ref = 300.0_R8   ! used for rho only; T stored per cell

    p0    = 1.0e5_R8
    R_air = 287.0_R8
    rho   = p0 / (R_air * T_ref)
    cdth  = cos(dth * 0.5_R8)
    sdth  = sin(dth * 0.5_R8)

    allocate(IOgas%block(1))
    IOgas%block(1)%Ni   = Ni
    IOgas%block(1)%Nj   = Nj
    IOgas%block(1)%Nk   = 1
    IOgas%block(1)%name = 'Block1'

    allocate(IOgas%block(1)%mesh(1:3, 0:Ni, 0:Nj, 0:1))
    allocate(IOgas%block(1)%vars(1:8, 1:Ni, 1:Nj, 1:1))

    !! Wedge mesh nodes
    do k = 0, 1; do j = 0, Nj; do i = 0, Ni
      rr = max(real(j,R8) * R0 / real(Nj,R8), r_min)
      IOgas%block(1)%mesh(1, i, j, k) = real(i,R8) * L0 / real(Ni,R8)
      IOgas%block(1)%mesh(2, i, j, k) = rr * cdth
      IOgas%block(1)%mesh(3, i, j, k) = rr * sdth * real(2*k - 1, R8)
    end do; end do; end do

    !! Cell-centre fields: non-uniform T from Tfun
    do k = 1, 1; do j = 1, Nj; do i = 1, Ni
      x_c = (real(i,R8) - 0.5_R8) * L0 / real(Ni,R8)
      r_c = (real(j,R8) - 0.5_R8) * R0 / real(Nj,R8)
      IOgas%block(1)%vars(1, i, j, k) = rho
      IOgas%block(1)%vars(2, i, j, k) = 0.0_R8
      IOgas%block(1)%vars(3, i, j, k) = 0.0_R8
      IOgas%block(1)%vars(4, i, j, k) = 0.0_R8
      IOgas%block(1)%vars(5, i, j, k) = p0
      IOgas%block(1)%vars(6, i, j, k) = Tfun(x_c, r_c)
      IOgas%block(1)%vars(7, i, j, k) = 1.4_R8
      IOgas%block(1)%vars(8, i, j, k) = R_air
    end do; end do; end do
  end subroutine Create_IOgas


  !! ----------------------------------------------------------
  !! Set boundary conditions with face-specific wall temperatures.
  subroutine Setup_BC(grid)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    integer :: b, ifa

    do b = 1, grid%nb
      do ifa = 1, 6
        grid%blk(b)%face(ifa)%Q   (:,:,:) = 0.0_R8
        grid%blk(b)%face(ifa)%Qtot(:,:)   = 0.0_R8
        grid%blk(b)%face(ifa)%a   (:,:,:) = 1.0_R8
      end do

      !! Endcap at x=0
      grid%blk(b)%face(1)%bc (:,:) = 301
      grid%blk(b)%face(1)%T  (:,:) = Tw_face1
      grid%blk(b)%face(1)%eps(:,:) = 1.0_R8

      !! Endcap at x=L
      grid%blk(b)%face(2)%bc (:,:) = 301
      grid%blk(b)%face(2)%T  (:,:) = Tw_face2
      grid%blk(b)%face(2)%eps(:,:) = 1.0_R8

      !! Axis (r=0) — specular reflection
      grid%blk(b)%face(3)%bc (:,:) = 300
      grid%blk(b)%face(3)%T  (:,:) = 0.0_R8
      grid%blk(b)%face(3)%eps(:,:) = 0.0_R8

      !! Lateral wall (r=R)
      grid%blk(b)%face(4)%bc (:,:) = 301
      grid%blk(b)%face(4)%T  (:,:) = Tw_lat
      grid%blk(b)%face(4)%eps(:,:) = 1.0_R8

      !! Azimuthal faces — periodic/axisymmetric
      grid%blk(b)%face(5)%bc (:,:) = 200
      grid%blk(b)%face(5)%T  (:,:) = 0.0_R8
      grid%blk(b)%face(5)%eps(:,:) = 0.0_R8

      grid%blk(b)%face(6)%bc (:,:) = 200
      grid%blk(b)%face(6)%T  (:,:) = 0.0_R8
      grid%blk(b)%face(6)%eps(:,:) = 0.0_R8
    end do
  end subroutine Setup_BC


  !! ----------------------------------------------------------
  !! Set radiative properties with non-uniform Ib(i,j) from Tfun.
  subroutine Setup_Rad_NonUniform(grid, Ni, Nj, R0, L0)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    integer,  intent(in) :: Ni, Nj
    real(R8), intent(in) :: R0, L0
    integer  :: b, i, j, k, Nx, Ny, Nz
    real(R8) :: x_c, r_c, T_c

    do b = 1, grid%nb
      Nx = grid%blk(b)%dim(1)
      Ny = grid%blk(b)%dim(2)
      Nz = grid%blk(b)%dim(3)

      allocate(grid%blk(b)%ka    (1:Nx, 1:Ny, 1:Nz, 0:Ngg))
      allocate(grid%blk(b)%a     (1:Nx, 1:Ny, 1:Nz, 0:Ngg))
      allocate(grid%blk(b)%Ib    (1:Nx, 1:Ny, 1:Nz))
      allocate(grid%blk(b)%source(1:Nx, 1:Ny, 1:Nz))
      allocate(grid%blk(b)%G     (1:Nx, 1:Ny, 1:Nz))

      grid%blk(b)%ka    (:,:,:,:) = k_user
      grid%blk(b)%a     (:,:,:,:) = 1.0_R8
      grid%blk(b)%source(:,:,:)   = 0.0_R8
      grid%blk(b)%G     (:,:,:)   = 0.0_R8

      do k = 1, Nz; do j = 1, Ny; do i = 1, Nx
        x_c = (real(i,R8) - 0.5_R8) * L0 / real(Ni,R8)
        r_c = (real(j,R8) - 0.5_R8) * R0 / real(Nj,R8)
        T_c = Tfun(x_c, r_c)
        grid%blk(b)%Ib(i, j, k) = sigma * T_c**4 / pi
      end do; end do; end do
    end do
  end subroutine Setup_Rad_NonUniform


  !! ----------------------------------------------------------
  !! Write source term at two radial positions and lateral wall flux vs x.
  subroutine Write_Output(grid, Ni, L0)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(in) :: grid
    integer,  intent(in) :: Ni
    real(R8), intent(in) :: L0
    integer :: b, i, unit1, unit2, unit3
    real(R8) :: dx, x

    b  = 1
    dx = L0 / real(Ni, R8)

    call execute_command_line('mkdir -p OUTPUT', wait=.true.)
    open(newunit=unit1, file='OUTPUT/source_j25.dat')
    write(unit1,'(A,F6.4,A,F5.1,A)') &
      '# x [m]              source [W/m3]   r_c=', &
      (real(j_25,R8)-0.5_R8)*R_cyl/real(Nj_rad,R8), ' m (', &
      (real(j_25,R8)-0.5_R8)/real(Nj_rad,R8)*100.0_R8, '% R)'
    do i = 1, Ni
      x = (real(i,R8) - 0.5_R8) * dx
      write(unit1,'(2ES20.8)') x, grid%blk(b)%source(i, j_25, 1)
    end do
    close(unit1)

    open(newunit=unit2, file='OUTPUT/source_j75.dat')
    write(unit2,'(A,F6.4,A,F5.1,A)') &
      '# x [m]              source [W/m3]   r_c=', &
      (real(j_75,R8)-0.5_R8)*R_cyl/real(Nj_rad,R8), ' m (', &
      (real(j_75,R8)-0.5_R8)/real(Nj_rad,R8)*100.0_R8, '% R)'
    do i = 1, Ni
      x = (real(i,R8) - 0.5_R8) * dx
      write(unit2,'(2ES20.8)') x, grid%blk(b)%source(i, j_75, 1)
    end do
    close(unit2)

    open(newunit=unit3, file='OUTPUT/heatflux_lat.dat')
    write(unit3,'(A)') '# x [m]              Qtot [W/m2]     (lateral wall, r=R)'
    do i = 1, Ni
      x = (real(i,R8) - 0.5_R8) * dx
      write(unit3,'(2ES20.8)') x, grid%blk(b)%face(4)%Qtot(i, 1)
    end do
    close(unit3)
  end subroutine Write_Output

end program furnace_test
