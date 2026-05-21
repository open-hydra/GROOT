!! ============================================================
!! Homogeneous cube radiative transfer test.
!!
!! Geometry : 1 x 1 x 1 m³ cube, 30x30x30 cells.
!! Medium   : uniform T = 1000 K, absorption coeff k = k_user.
!! Walls    : black (eps=1), cold (Tw = 0 K).
!! Model    : DTM, spectral model = const.
!!
!! The test is self-contained — no external file I/O is needed
!! for setup (mesh and field are built in memory).
!! After solving, the z-axis centerline source term and wall heat
!! flux on face 1 are written to OUTPUT/ with the k value in the
!! file name (e.g. source_centerline_k1.dat) for comparison with
!! reference solutions in data/.
!! ============================================================
program homcube_test
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

  type(GROOT_simulation_type) :: simulation
  type(orion_data)            :: IOgas
  real(R8) :: t1, t2
  integer  :: b, i, j, k, ifa, i1, i2, n1, n2

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
  !! Configure via typed config objects (no input.ini needed)
  obj_dtm%Nr             = 256
  obj_dtm%tol            = 1.0e-3_R8
  obj_dtm%itermax        = 1
  obj_dtm%source         = .true.
  obj_dtm%twoDax         = .false.
  obj_dtm%wedge_angle_deg= 1.0_R8
  obj_dtm%optically_thin = .false.

  obj_rad_model%model    = 'const'
  obj_rad_model%eps_wall = 1.0_R8
  obj_rad_model%k_user   = 1.0_R8

  obj_io%sol_format = 'tecplot ascii'
  obj_io%gas_format = 'tecplot ascii'


  !! Set species manually (N2 inert, nscrad=0 for const model)
  nsc  = 1
  nsoot= 0
  nrans = 0
  allocate(cfd_species(1:nsc), wm_tab(1:nsc))
  cfd_species(1) = 'N2'
  wm_tab(1)      = 28.0_R8
  nscrad = 0
  nscdat = 0
  allocate(rad_species(0), ind_species(0))

  !! Sync to GROOT_Mod_Parameters globals (needed by physics modules)
  call Assign_Setup()

  write(*,'(A)')
  write(*,'(A)') ' ========================================'
  write(*,'(A)') ' Homogeneous cube test'
  write(*,'(A,I0)') '   Nr    = ', Nr
  write(*,'(A,F6.2)') '   k_user = ', k_user
  write(*,'(A)') ' ========================================'
  write(*,'(A)')

  !! --------------------------------------------------------
  !! Build IOgas in memory: 30x30x30 cube, 1x1x1 m³, T=1000 K
  call Create_IOgas(IOgas, nx=30, ny=30, nz=30, &
                    Lx=1.0_R8, Ly=1.0_R8, Lz=1.0_R8, &
                    T_gas=1000.0_R8)

  !! --------------------------------------------------------
  !! Allocate grid data structures and compute mesh metrics
  call Setup_Data_Structure(simulation%domain, IOgas)
  call Setup_Metrics(simulation%domain, IOgas)

  !! --------------------------------------------------------
  !! Set wall BCs: bc=6 (opaque wall), T=0, eps=1 on all faces
  call Setup_Wall_Direct(simulation%domain, Tw=0.0_R8, eps=1.0_R8)

  !! --------------------------------------------------------
  !! Set radiative properties for const model
  call Setup_Rad_Const(simulation%domain, T_gas=1000.0_R8)

  !! --------------------------------------------------------
  !! Solve
  call cpu_time(t1)
  call GROOT_solve(simulation)
  call cpu_time(t2)

  write(*,'(A)')
  write(*,'(A,F8.4,A)') '   CPU time: ', &
    (t2-t1)/real(obj_sim_param%nthreads,R8)/60.0_R8, ' min'

  !! --------------------------------------------------------
  !! Write output for post-processing
  call Write_Output(simulation%domain)

  write(*,'(A)') ' Done. See OUTPUT/source_centerline_k*.dat and OUTPUT/heatflux_face1_k*.dat'
  write(*,'(A)')

contains

  !! ----------------------------------------------------------
  subroutine Create_IOgas(IOgas, nx, ny, nz, Lx, Ly, Lz, T_gas)
    implicit none
    type(orion_data), intent(out) :: IOgas
    integer,  intent(in) :: nx, ny, nz
    real(R8), intent(in) :: Lx, Ly, Lz, T_gas
    integer :: i, j, k
    real(R8) :: p0, R_air, rho

    p0    = 1.0e5_R8
    R_air = 287.0_R8
    rho   = p0 / (R_air * T_gas)

    allocate(IOgas%block(1))
    IOgas%block(1)%Ni = nx
    IOgas%block(1)%Nj = ny
    IOgas%block(1)%Nk = nz
    IOgas%block(1)%name = 'Block1'

    allocate(IOgas%block(1)%mesh(1:3, 0:nx, 0:ny, 0:nz))
    allocate(IOgas%block(1)%vars(1:8, 1:nx, 1:ny, 1:nz))

    !! Uniform Cartesian mesh
    do k = 0, nz; do j = 0, ny; do i = 0, nx
      IOgas%block(1)%mesh(1,i,j,k) = real(i,R8) * Lx / real(nx,R8)
      IOgas%block(1)%mesh(2,i,j,k) = real(j,R8) * Ly / real(ny,R8)
      IOgas%block(1)%mesh(3,i,j,k) = real(k,R8) * Lz / real(nz,R8)
    end do; end do; end do

    !! Uniform field: [rho, u, v, w, p, T, gamma, R]
    !! with nsc=1, nrans=0, nsoot=0:
    !!   vars(nsc)      = rho(N2)       index 1
    !!   vars(nsc+4)    = p             index 5
    !!   vars(nsc+5)    = T             index 6
    !!   vars(nsc+7)    = R             index 8
    do k = 1, nz; do j = 1, ny; do i = 1, nx
      IOgas%block(1)%vars(1, i,j,k) = rho       ! rho_N2
      IOgas%block(1)%vars(2, i,j,k) = 0.0_R8    ! u
      IOgas%block(1)%vars(3, i,j,k) = 0.0_R8    ! v
      IOgas%block(1)%vars(4, i,j,k) = 0.0_R8    ! w
      IOgas%block(1)%vars(5, i,j,k) = p0        ! p
      IOgas%block(1)%vars(6, i,j,k) = T_gas     ! T
      IOgas%block(1)%vars(7, i,j,k) = 1.4_R8    ! gamma
      IOgas%block(1)%vars(8, i,j,k) = R_air     ! R
    end do; end do; end do
  end subroutine Create_IOgas


  !! ----------------------------------------------------------
  !! Directly set BCs on all faces: bc=bcval, T=Tw, eps=eps_val
  subroutine Setup_Wall_Direct(grid, Tw, eps)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    real(R8), intent(in) :: Tw, eps
    integer :: b, ifa

    do b = 1, grid%nb
      do ifa = 1, 6
        grid%blk(b)%face(ifa)%bc  (:,:)   = 301
        grid%blk(b)%face(ifa)%T  (:,:)   = Tw
        grid%blk(b)%face(ifa)%eps(:,:)   = eps
        grid%blk(b)%face(ifa)%Q  (:,:,:) = 0.0_R8
        grid%blk(b)%face(ifa)%Qtot(:,:)  = 0.0_R8
        grid%blk(b)%face(ifa)%a  (:,:,:) = 1.0_R8
      end do
    end do
  end subroutine Setup_Wall_Direct


  !! ----------------------------------------------------------
  !! Set radiative properties for const model (no species data needed)
  subroutine Setup_Rad_Const(grid, T_gas)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    real(R8), intent(in) :: T_gas
    integer :: b, Nx, Ny, Nz
    real(R8) :: Ib_val

    Ib_val = sigma * T_gas**4 / pi   ! Blackbody intensity [W/m²/sr]

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
      grid%blk(b)%Ib    (:,:,:)   = Ib_val
      grid%blk(b)%source(:,:,:)   = 0.0_R8
      grid%blk(b)%G     (:,:,:)   = 0.0_R8
    end do
  end subroutine Setup_Rad_Const


  !! ----------------------------------------------------------
  !! Write z-axis centerline source term and face-1 wall heat flux.
  !! File names include the absorption coefficient, e.g. source_centerline_k1.dat
  subroutine Write_Output(grid)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(in) :: grid
    integer :: b, i, j, k, unit1, unit2
    integer :: Nx, Ny, Nz, ic, jc
    real(R8) :: dz, z
    character(len=8) :: k_str

    !! Build k label: 0.1 → "0.1", 1.0 → "1", 10.0 → "10"
    if (abs(k_user - 0.1_R8) < 1.0e-6_R8) then
      k_str = '0.1'
    elseif (abs(k_user - 10.0_R8) < 1.0e-6_R8) then
      k_str = '10'
    else
      k_str = '1'
    end if

    b  = 1
    Nx = grid%blk(b)%dim(1)
    Ny = grid%blk(b)%dim(2)
    Nz = grid%blk(b)%dim(3)
    ic = Nx/2 + 1
    jc = Ny/2 + 1
    dz = grid%blk(b)%z(1,1,Nz+1) - grid%blk(b)%z(1,1,1)

    open(newunit=unit1, file='OUTPUT/source_centerline_k'//trim(k_str)//'.dat')
    write(unit1,'(A)') '# z   source [W/m³]'
    do k = 1, Nz
      z = (real(k,R8) - 0.5_R8) * dz / real(Nz,R8)
      write(unit1,'(2ES20.8)') z, grid%blk(b)%source(ic, jc, k)
    end do
    close(unit1)

    open(newunit=unit2, file='OUTPUT/heatflux_face1_k'//trim(k_str)//'.dat')
    write(unit2,'(A)') '# z   Qtot [W/m²]  (face 1, j=jc)'
    do k = 1, Nz
      z = (real(k,R8) - 0.5_R8) * dz / real(Nz,R8)
      write(unit2,'(2ES20.8)') z, grid%blk(b)%face(1)%Qtot(jc, k)
    end do
    close(unit2)
  end subroutine Write_Output

end program homcube_test
