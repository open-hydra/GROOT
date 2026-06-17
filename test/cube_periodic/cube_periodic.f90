!! ============================================================
!! Symmetry (periodic) BC test for the homogeneous cube.
!!
!! Geometry : 0.5 x 0.5 x 0.5 m³ cube (1/8 of the reference cube),
!!            15 x 15 x 15 cells (same dx = 1/30 m as homcube).
!! Medium   : uniform T = 1000 K, absorption coeff k = k_user.
!! BCs      :
!!   Faces 1, 3, 5  (x=0, y=0, z=0) : opaque black wall  (bc=302, Tw=0, eps=1)
!!   Faces 2, 4, 6  (x=0.5, y=0.5, z=0.5) : symmetry plane (bc=300)
!!
!! Expected result:
!!   By symmetry of the full 1x1x1 cube, the radiation field in [0,0.5]^3
!!   is identical to the corresponding octant of the full cube.
!!   The source term and wall flux must therefore match the homcube reference
!!   within the discretisation error.  Since cell spacing is identical
!!   (dx=1/30 m), cell (i,j,k) for i,j,k<=15 maps one-to-one.
!!
!!   Comparison:  source at (ic=Nx, jc=Ny, k=1..Nz)  [near-symmetry corner,
!!   x=y≈0.483 m]  vs  ref_source_k1_T_1000.dat  [computed at x=y=0.5 m].
!!   The two positions are mirror images about x=0.5, so by symmetry they
!!   are equal to within the angular-quadrature error of the reference.
!! ============================================================
program cube_periodic_test
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
  !! Configure (same physics as homcube)
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


  !! Set species manually
  nsc  = 1
  nrans   = 0
  nsoot= 0
  allocate(cfd_species(1:nsc), wm_tab(1:nsc))
  cfd_species(1) = 'N2'
  wm_tab(1)      = 28.0_R8
  nscrad = 0
  nscdat = 0
  allocate(rad_species(0), ind_species(0))

  call Assign_Setup()

  write(*,'(A)')
  write(*,'(A)') ' =============================================='
  write(*,'(A)') ' Cube periodic (symmetry BC) test'
  write(*,'(A)') '   Domain  : [0, 0.5]^3 m   (1/8 of full cube)'
  write(*,'(A)') '   Mesh    : 15 x 15 x 15   (dx = 1/30 m)'
  write(*,'(A)') '   bc 1,3,5: black wall (bc=302, Tw=0, eps=1)'
  write(*,'(A)') '   bc 2,4,6: symmetry plane (bc=300)'
  write(*,'(A,I0)') '   Nr      = ', Nr
  write(*,'(A,F6.2)') '   k_user  = ', k_user
  write(*,'(A)') ' =============================================='
  write(*,'(A)')

  !! --------------------------------------------------------
  !! Build IOgas in memory: 15x15x15, 0.5x0.5x0.5 m, T=1000 K
  call Create_IOgas(IOgas, nx=15, ny=15, nz=15, &
                    Lx=0.5_R8, Ly=0.5_R8, Lz=0.5_R8, &
                    T_gas=1000.0_R8)

  !! --------------------------------------------------------
  !! Allocate grid and compute mesh metrics
  call Setup_Data_Structure(simulation%domain, IOgas)
  call Setup_Metrics(simulation%domain, IOgas)

  !! --------------------------------------------------------
  !! Set BCs:
  !!   odd faces (1,3,5) → black wall bc=302, T=0, eps=1
  !!   even faces (2,4,6) → symmetry bc=300
  call Setup_BC(simulation%domain, Tw=0.0_R8, eps=1.0_R8)

  !! --------------------------------------------------------
  !! Set radiative properties (const model)
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
  call Write_Output(simulation%domain)

  write(*,'(A)') ' Done. See OUTPUT/source_corner.dat and OUTPUT/heatflux_face1.dat'
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

    do k = 0, nz; do j = 0, ny; do i = 0, nx
      IOgas%block(1)%mesh(1,i,j,k) = real(i,R8) * Lx / real(nx,R8)
      IOgas%block(1)%mesh(2,i,j,k) = real(j,R8) * Ly / real(ny,R8)
      IOgas%block(1)%mesh(3,i,j,k) = real(k,R8) * Lz / real(nz,R8)
    end do; end do; end do

    do k = 1, nz; do j = 1, ny; do i = 1, nx
      IOgas%block(1)%vars(1, i,j,k) = rho
      IOgas%block(1)%vars(2, i,j,k) = 0.0_R8
      IOgas%block(1)%vars(3, i,j,k) = 0.0_R8
      IOgas%block(1)%vars(4, i,j,k) = 0.0_R8
      IOgas%block(1)%vars(5, i,j,k) = p0
      IOgas%block(1)%vars(6, i,j,k) = T_gas
      IOgas%block(1)%vars(7, i,j,k) = 1.4_R8
      IOgas%block(1)%vars(8, i,j,k) = R_air
    end do; end do; end do
  end subroutine Create_IOgas


  !! ----------------------------------------------------------
  !! Set BCs:
  !!   faces 1, 3, 5 (x=0, y=0, z=0) → opaque black wall (bc=302)
  !!   faces 2, 4, 6 (x=Lx, y=Ly, z=Lz) → symmetry (bc=300)
  subroutine Setup_BC(grid, Tw, eps)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    real(R8), intent(in) :: Tw, eps
    integer :: b, ifa

    do b = 1, grid%nb
      do ifa = 1, 6
        !! Initialise flux and weight arrays on all faces
        grid%blk(b)%face(ifa)%Q   (:,:,:) = 0.0_R8
        grid%blk(b)%face(ifa)%Qtot(:,:)   = 0.0_R8
        grid%blk(b)%face(ifa)%a   (:,:,:) = 1.0_R8

        if (mod(ifa, 2) == 1) then
          !! Odd face → black cold wall
          grid%blk(b)%face(ifa)%bc (:,:) = 301
          grid%blk(b)%face(ifa)%T  (:,:) = Tw
          grid%blk(b)%face(ifa)%eps(:,:) = eps
        else
          !! Even face → symmetry plane
          grid%blk(b)%face(ifa)%bc (:,:) = 300
          grid%blk(b)%face(ifa)%T  (:,:) = 0.0_R8
          grid%blk(b)%face(ifa)%eps(:,:) = 0.0_R8
        end if
      end do
    end do
  end subroutine Setup_BC


  !! ----------------------------------------------------------
  subroutine Setup_Rad_Const(grid, T_gas)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    real(R8), intent(in) :: T_gas
    integer :: b, Nx, Ny, Nz
    real(R8) :: Ib_val

    Ib_val = sigma * T_gas**4 / pi

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
  !! Write source at the near-symmetry corner (ic=Nx, jc=Ny) vs z,
  !! and wall heat flux at face 1 (x=0) at jc=Ny vs z.
  !!
  !! These positions (x=y≈0.483 m) are mirror images of the homcube
  !! centerline (x=y≈0.517 m) about the symmetry planes, so by
  !! problem symmetry the values must match the homcube reference.
  subroutine Write_Output(grid)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(in) :: grid
    integer :: b, k, unit1, unit2
    integer :: Nx, Ny, Nz, ic, jc
    real(R8) :: dz, z

    b  = 1
    Nx = grid%blk(b)%dim(1)
    Ny = grid%blk(b)%dim(2)
    Nz = grid%blk(b)%dim(3)
    ic = Nx           ! near symmetry corner: i = Nx  (x ≈ 0.483 m)
    jc = Ny           !                       j = Ny  (y ≈ 0.483 m)
    dz = grid%blk(b)%z(1,1,Nz+1) - grid%blk(b)%z(1,1,1)

    call execute_command_line('mkdir -p OUTPUT', wait=.true.)
    open(newunit=unit1, file='OUTPUT/source_corner.dat')
    write(unit1,'(A)') '# z [m]              source [W/m3]'
    write(unit1,'(A)') '# position: x=y≈0.483 m (near symmetry planes, same as'
    write(unit1,'(A)') '# mirror of homcube centerline at x=y≈0.517 m)'
    do k = 1, Nz
      z = (real(k,R8) - 0.5_R8) * dz / real(Nz,R8)
      write(unit1,'(2ES20.8)') z, grid%blk(b)%source(ic, jc, k)
    end do
    close(unit1)

    open(newunit=unit2, file='OUTPUT/heatflux_face1.dat')
    write(unit2,'(A)') '# z [m]              Qtot [W/m2]  (face 1, x=0, jc=Ny)'
    do k = 1, Nz
      z = (real(k,R8) - 0.5_R8) * dz / real(Nz,R8)
      write(unit2,'(2ES20.8)') z, grid%blk(b)%face(1)%Qtot(jc, k)
    end do
    close(unit2)
  end subroutine Write_Output

end program cube_periodic_test
