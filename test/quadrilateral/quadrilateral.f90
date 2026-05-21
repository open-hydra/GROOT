!! ============================================================
!! Quadrilateral (scalene trapezoid) radiative transfer test.
!!
!! Geometry : scalene trapezoid ABCD (counter-clockwise),
!!            extruded as a thin slab of thickness Lz=0.01 m.
!!   A(0.0, 0.0)  — bottom-left
!!   B(2.2, 0.0)  — bottom-right
!!   C(1.5, 1.2)  — top-right
!!   D(0.5, 1.0)  — top-left
!!   Mesh: Ni×Nj×1 cells (bilinear transfinite interpolation).
!!
!! Medium   : uniform T=3000 K, absorption coeff k=k_user.
!! BCs      :
!!   Faces 1–4 (quadrilateral walls) → black walls (bc=301, Tw=0)
!!   Faces 5,6 (z-slab faces)        → specular   (bc=300)
!!
!! Output   : wall heat flux on face 3 (bottom AB) vs x  [W/m²]
!!
!! Reference: data/ref_k{k}.dat (from old GROOT; note the old
!!            code segfaulted on this case, so ref is indicative).
!! ============================================================
program quadrilateral_test
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

  !! Mesh parameters
  integer,  parameter :: Ni = 19          !< cells in i (horizontal)
  integer,  parameter :: Nj = 19          !< cells in j (vertical)
  real(R8), parameter :: Lz = 0.01_R8    !< slab thickness [m]

  !! Quadrilateral vertices (counter-clockwise: A→B→C→D)
  real(R8), parameter :: Ax=0.0_R8, Ay=0.0_R8   !< A: bottom-left
  real(R8), parameter :: Bx=2.2_R8, By=0.0_R8   !< B: bottom-right
  real(R8), parameter :: Cx=1.5_R8, Cy=1.2_R8   !< C: top-right
  real(R8), parameter :: Dx=0.5_R8, Dy=1.0_R8   !< D: top-left

  real(R8), parameter :: T_gas = 3000.0_R8       !< gas temperature [K]

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
  !! Configure (3-D slab — NOT axisymmetric)
  obj_dtm%Nr             = 256
  obj_dtm%tol            = 1.0e-3_R8
  obj_dtm%itermax        = 1
  obj_dtm%source         = .false.
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
  write(*,'(A)') ' Quadrilateral (scalene trapezoid) test'
  write(*,'(A)') '   A(0.0, 0.0)  B(2.2, 0.0)'
  write(*,'(A)') '   D(0.5, 1.0)  C(1.5, 1.2)'
  write(*,'(A,I0,A,I0)') '   Mesh   : ', Ni, ' x ', Nj, ' x 1  (bilinear TFI)'
  write(*,'(A,F5.3,A)') '   Lz     = ', Lz, ' m (thin slab)'
  write(*,'(A,F8.1,A)') '   T_gas  = ', T_gas, ' K'
  write(*,'(A,F6.2)') '   k_user = ', k_user
  write(*,'(A,I0)') '   Nr     = ', Nr
  write(*,'(A)') ' =============================================='
  write(*,'(A)')

  !! --------------------------------------------------------
  !! Build IOgas in memory: bilinear TFI mesh, uniform T
  call Create_IOgas(IOgas, Ni, Nj, T_gas)

  !! --------------------------------------------------------
  !! Allocate grid and compute mesh metrics
  call Setup_Data_Structure(simulation%domain, IOgas)
  call Setup_Metrics(simulation%domain, IOgas)

  !! --------------------------------------------------------
  !! Set BCs
  call Setup_BC(simulation%domain)

  !! --------------------------------------------------------
  !! Set radiative properties
  call Setup_Rad_Const(simulation%domain, T_gas)

  !! --------------------------------------------------------
  !! Solve
  call cpu_time(t1)
  call GROOT_solve(simulation)
  call cpu_time(t2)

  write(*,'(A)')
  write(*,'(A,F8.4,A)') '   CPU time: ', &
    (t2-t1)/real(obj_sim_param%nthreads,R8)/60.0_R8, ' min'

  !! --------------------------------------------------------
  call Write_Output(simulation%domain, Ni)

  write(*,'(A)') ' Done. See OUTPUT/heatflux_face3.dat'
  write(*,'(A)')

contains

  !! ----------------------------------------------------------
  !! Build IOgas with bilinear TFI mesh for the quadrilateral.
  !!
  !! Parametric mapping (xi in [0,1], eta in [0,1]):
  !!   P(xi,eta) = (1-xi)(1-eta)*A + xi(1-eta)*B + xi*eta*C + (1-xi)*eta*D
  !!
  !! Node (i,j): xi = i/Ni,  eta = j/Nj  (i=0..Ni, j=0..Nj)
  !! z nodes:    k=0 → -Lz,  k=1 → +Lz
  subroutine Create_IOgas(IOgas, Nii, Njj, T0)
    implicit none
    type(orion_data), intent(out) :: IOgas
    integer,  intent(in) :: Nii, Njj
    real(R8), intent(in) :: T0
    integer  :: i, j, k
    real(R8) :: p0, R_air, rho, xi, eta

    p0    = 1.0e5_R8
    R_air = 287.0_R8
    rho   = p0 / (R_air * T0)

    allocate(IOgas%block(1))
    IOgas%block(1)%Ni   = Nii
    IOgas%block(1)%Nj   = Njj
    IOgas%block(1)%Nk   = 1
    IOgas%block(1)%name = 'Block1'

    allocate(IOgas%block(1)%mesh(1:3, 0:Nii, 0:Njj, 0:1))
    allocate(IOgas%block(1)%vars(1:8, 1:Nii, 1:Njj, 1:1))

    !! Bilinear TFI mesh nodes
    do k = 0, 1; do j = 0, Njj; do i = 0, Nii
      xi  = real(i, R8) / real(Nii, R8)
      eta = real(j, R8) / real(Njj, R8)
      IOgas%block(1)%mesh(1, i, j, k) = (1.0_R8-xi)*(1.0_R8-eta)*Ax &
                                        +            xi*(1.0_R8-eta)*Bx &
                                        +                     xi*eta*Cx &
                                        +       (1.0_R8-xi)*eta*Dx
      IOgas%block(1)%mesh(2, i, j, k) = (1.0_R8-xi)*(1.0_R8-eta)*Ay &
                                        +            xi*(1.0_R8-eta)*By &
                                        +                     xi*eta*Cy &
                                        +       (1.0_R8-xi)*eta*Dy
      IOgas%block(1)%mesh(3, i, j, k) = Lz * real(2*k - 1, R8)
    end do; end do; end do

    !! Uniform field: T = T_gas
    do k = 1, 1; do j = 1, Njj; do i = 1, Nii
      IOgas%block(1)%vars(1, i, j, k) = rho
      IOgas%block(1)%vars(2, i, j, k) = 0.0_R8
      IOgas%block(1)%vars(3, i, j, k) = 0.0_R8
      IOgas%block(1)%vars(4, i, j, k) = 0.0_R8
      IOgas%block(1)%vars(5, i, j, k) = p0
      IOgas%block(1)%vars(6, i, j, k) = T0
      IOgas%block(1)%vars(7, i, j, k) = 1.4_R8
      IOgas%block(1)%vars(8, i, j, k) = R_air
    end do; end do; end do
  end subroutine Create_IOgas


  !! ----------------------------------------------------------
  !! Set BCs: faces 1–4 are opaque black walls (Tw=0),
  !!          faces 5,6 are specular (z-symmetry of thin slab).
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

      !! Quadrilateral side walls: black, cold
      do ifa = 1, 4
        grid%blk(b)%face(ifa)%bc (:,:) = 301
        grid%blk(b)%face(ifa)%T  (:,:) = 0.0_R8
        grid%blk(b)%face(ifa)%eps(:,:) = 1.0_R8
      end do

      !! z-slab faces: specular (problem is 2D by symmetry)
      grid%blk(b)%face(5)%bc (:,:) = 300
      grid%blk(b)%face(5)%T  (:,:) = 0.0_R8
      grid%blk(b)%face(5)%eps(:,:) = 0.0_R8

      grid%blk(b)%face(6)%bc (:,:) = 300
      grid%blk(b)%face(6)%T  (:,:) = 0.0_R8
      grid%blk(b)%face(6)%eps(:,:) = 0.0_R8
    end do
  end subroutine Setup_BC


  !! ----------------------------------------------------------
  subroutine Setup_Rad_Const(grid, T0)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    real(R8), intent(in) :: T0
    integer :: b, Nx, Ny, Nz
    real(R8) :: Ib_val

    Ib_val = sigma * T0**4 / pi

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
  !! Write wall heat flux on face 3 (bottom wall AB, j=1) vs
  !! x-coordinate of each cell centre along that face.
  !!
  !! Face 3 (j-min) spans i=1..Ni cells; Qtot(i,1) with Nk=1.
  !! Bottom edge is straight A→B, so x_c = Bx*(i-0.5)/Ni.
  subroutine Write_Output(grid, Nii)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(in) :: grid
    integer, intent(in) :: Nii
    integer :: b, i, unit1
    real(R8) :: x_c

    b = 1
    open(newunit=unit1, file='OUTPUT/heatflux_face3.dat')
    write(unit1,'(A)') '# x [m]              Qtot [W/m2]     (face 3, bottom AB)'
    do i = 1, Nii
      x_c = Bx * (real(i, R8) - 0.5_R8) / real(Nii, R8)
      write(unit1,'(2ES20.8)') x_c, grid%blk(b)%face(3)%Qtot(i, 1)
    end do
    close(unit1)
  end subroutine Write_Output

end program quadrilateral_test
