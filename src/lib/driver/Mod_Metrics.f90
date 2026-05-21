module GROOT_Mod_Metrics
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private
  public :: Setup_Metrics
  public :: Wedge
  public :: Import_Face_Coords
  public :: meshType

  integer  :: meshType   !< 1 = 1D, 2 = 2D, 3 = 3D


contains

  subroutine Setup_Metrics(grid, IOfield)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use Lib_ORION_Data
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    type(orion_data),        intent(inout) :: IOfield

    call Import_Nodes(grid, IOfield)
    call Import_Face_Coords(grid)
    call Check_Mesh_Type(grid)
  end subroutine Setup_Metrics


  subroutine Import_Nodes(grid, IOfield)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use Lib_ORION_Data
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    type(orion_data),        intent(inout) :: IOfield
    integer :: b

    grid%nb = size(IOfield%block)

    do b = 1, grid%nb
      grid%blk(b)%dim(1) = IOfield%block(b)%Ni
      grid%blk(b)%dim(2) = IOfield%block(b)%Nj
      grid%blk(b)%dim(3) = IOfield%block(b)%Nk

      associate(ni => grid%blk(b)%dim(1), &
                nj => grid%blk(b)%dim(2), &
                nk => grid%blk(b)%dim(3))
        allocate(grid%blk(b)%x(1:ni+1, 1:nj+1, 1:nk+1))
        allocate(grid%blk(b)%y(1:ni+1, 1:nj+1, 1:nk+1))
        allocate(grid%blk(b)%z(1:ni+1, 1:nj+1, 1:nk+1))

        grid%blk(b)%x = IOfield%block(b)%mesh(1, 0:ni, 0:nj, 0:nk)
        grid%blk(b)%y = IOfield%block(b)%mesh(2, 0:ni, 0:nj, 0:nk)
        grid%blk(b)%z = IOfield%block(b)%mesh(3, 0:ni, 0:nj, 0:nk)
      end associate
    end do
  end subroutine Import_Nodes


  !! @brief Copy nodal coordinates from block x/y/z arrays into face coordinate arrays.
  !! Must be called after Import_Nodes (which fills grid%blk%x/y/z) and after
  !! Setup_Data_Structure (which allocates face%x/y/z via face%alloc).
  subroutine Import_Face_Coords(grid)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    integer :: b, Nx, Ny, Nz

    do b = 1, grid%nb
      Nx = grid%blk(b)%dim(1)
      Ny = grid%blk(b)%dim(2)
      Nz = grid%blk(b)%dim(3)

      !! Face 1: i-min  (n=[Ny,Nz], node index i=1)
      grid%blk(b)%face(1)%x(1:Ny+1, 1:Nz+1) = grid%blk(b)%x(1,    1:Ny+1, 1:Nz+1)
      grid%blk(b)%face(1)%y(1:Ny+1, 1:Nz+1) = grid%blk(b)%y(1,    1:Ny+1, 1:Nz+1)
      grid%blk(b)%face(1)%z(1:Ny+1, 1:Nz+1) = grid%blk(b)%z(1,    1:Ny+1, 1:Nz+1)

      !! Face 2: i-max  (node index i=Nx+1)
      grid%blk(b)%face(2)%x(1:Ny+1, 1:Nz+1) = grid%blk(b)%x(Nx+1, 1:Ny+1, 1:Nz+1)
      grid%blk(b)%face(2)%y(1:Ny+1, 1:Nz+1) = grid%blk(b)%y(Nx+1, 1:Ny+1, 1:Nz+1)
      grid%blk(b)%face(2)%z(1:Ny+1, 1:Nz+1) = grid%blk(b)%z(Nx+1, 1:Ny+1, 1:Nz+1)

      !! Face 3: j-min  (n=[Nx,Nz], node index j=1)
      grid%blk(b)%face(3)%x(1:Nx+1, 1:Nz+1) = grid%blk(b)%x(1:Nx+1, 1,    1:Nz+1)
      grid%blk(b)%face(3)%y(1:Nx+1, 1:Nz+1) = grid%blk(b)%y(1:Nx+1, 1,    1:Nz+1)
      grid%blk(b)%face(3)%z(1:Nx+1, 1:Nz+1) = grid%blk(b)%z(1:Nx+1, 1,    1:Nz+1)

      !! Face 4: j-max  (node index j=Ny+1)
      grid%blk(b)%face(4)%x(1:Nx+1, 1:Nz+1) = grid%blk(b)%x(1:Nx+1, Ny+1, 1:Nz+1)
      grid%blk(b)%face(4)%y(1:Nx+1, 1:Nz+1) = grid%blk(b)%y(1:Nx+1, Ny+1, 1:Nz+1)
      grid%blk(b)%face(4)%z(1:Nx+1, 1:Nz+1) = grid%blk(b)%z(1:Nx+1, Ny+1, 1:Nz+1)

      !! Face 5: k-min  (n=[Nx,Ny], node index k=1)
      grid%blk(b)%face(5)%x(1:Nx+1, 1:Ny+1) = grid%blk(b)%x(1:Nx+1, 1:Ny+1, 1   )
      grid%blk(b)%face(5)%y(1:Nx+1, 1:Ny+1) = grid%blk(b)%y(1:Nx+1, 1:Ny+1, 1   )
      grid%blk(b)%face(5)%z(1:Nx+1, 1:Ny+1) = grid%blk(b)%z(1:Nx+1, 1:Ny+1, 1   )

      !! Face 6: k-max  (node index k=Nz+1)
      grid%blk(b)%face(6)%x(1:Nx+1, 1:Ny+1) = grid%blk(b)%x(1:Nx+1, 1:Ny+1, Nz+1)
      grid%blk(b)%face(6)%y(1:Nx+1, 1:Ny+1) = grid%blk(b)%y(1:Nx+1, 1:Ny+1, Nz+1)
      grid%blk(b)%face(6)%z(1:Nx+1, 1:Ny+1) = grid%blk(b)%z(1:Nx+1, 1:Ny+1, Nz+1)
    end do
  end subroutine Import_Face_Coords


  subroutine Check_Mesh_Type(grid)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use GROOT_Global_m, only: wedge_angle_deg, pi
    implicit none
    type(GROOT_domain_type), intent(in) :: grid
    integer :: b, nk_total
    real(R8) :: delthe     !< Axisymmetric grid angle (rad)
    nk_total = 0
    do b = 1, grid%nb
      nk_total = nk_total + grid%blk(b)%dim(3)
    end do

    if (nk_total == grid%nb) then
      meshType = 2
      delthe   = wedge_angle_deg * pi / 180.0_R8
    else
      meshType = 3
      delthe   = 0.0_R8
    end if
  end subroutine Check_Mesh_Type


  !! @brief Apply 2-D axisymmetric wedge correction to node coordinates.
  !! Nodes at k=1 are placed at theta=-angle/2 and nodes at k=Nz+1 at theta=+angle/2,
  !! interpolating linearly between them.  The true radius is recovered as sqrt(y^2+z^2)
  !! so the routine is idempotent regardless of the angular offset stored in the input mesh.
  subroutine Wedge(grid, angle)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    real(R8),                intent(in)    :: angle
    integer  :: b, i, j, k
    real(R8) :: r, theta_k

    do b = 1, grid%nb
      associate(Nx => grid%blk(b)%dim(1), &
                Ny => grid%blk(b)%dim(2), &
                Nz => grid%blk(b)%dim(3))
        do k = 1, Nz+1; do j = 1, Ny+1; do i = 1, Nx+1
          r = sqrt(grid%blk(b)%y(i,j,k)**2 + grid%blk(b)%z(i,j,k)**2)
          theta_k = -angle/2.0_R8 + real(k-1, R8) * angle / real(Nz, R8)
          grid%blk(b)%y(i,j,k) =  r * cos(theta_k)
          grid%blk(b)%z(i,j,k) =  r * sin(theta_k)
        end do; end do; end do
      end associate
    end do
  end subroutine Wedge

end module GROOT_Mod_Metrics
