module GROOT_Mod_Allocate_Data
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private
  public :: Setup_Data_Structure
  public :: Deallocate_Data

contains

  subroutine Setup_Data_Structure(grid, IOfield)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use Lib_ORION_Data
    use GROOT_Mod_Face, only: SBface
    use GROOT_Global_m, only: nsc
    use GROOT_Mod_Phase, only: Setup_Face
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    type(orion_data),        intent(inout) :: IOfield
    integer :: ni, nj, nk, nbound
    integer :: b, ifa

    if (allocated(grid%blk)) deallocate(grid%blk)
    grid%nb = size(IOfield%block)
    allocate(grid%blk(1:grid%nb))

    nbound = 0
    do b = 1, grid%nb

      ni = IOfield%block(b)%Ni
      nj = IOfield%block(b)%Nj
      nk = IOfield%block(b)%Nk

      !! Face sizes
      nbound = nbound + 2*nj*nk + 2*ni*nk + 2*nj*ni

      grid%blk(b)%face(1)%n = [nj, nk]
      grid%blk(b)%face(2)%n = [nj, nk]
      grid%blk(b)%face(3)%n = [ni, nk]
      grid%blk(b)%face(4)%n = [ni, nk]
      grid%blk(b)%face(5)%n = [ni, nj]
      grid%blk(b)%face(6)%n = [ni, nj]


      do ifa = 1, 6
        call grid%blk(b)%face(ifa)%alloc
      end do

      call Setup_Face(b, 1, [nj, nk], grid)
      call Setup_Face(b, 2, [nj, nk], grid)
      call Setup_Face(b, 3, [ni, nk], grid)
      call Setup_Face(b, 4, [ni, nk], grid)
      call Setup_Face(b, 5, [ni, nj], grid)
      call Setup_Face(b, 6, [ni, nj], grid)

      !! Gaseous phase (gas_phase is a direct component, not allocatable)
      allocate(grid%blk(b)%gas_phase%rho(1:nsc,  1:ni, 1:nj, 1:nk))
      allocate(grid%blk(b)%gas_phase%T  (         1:ni, 1:nj, 1:nk))
      allocate(grid%blk(b)%gas_phase%p  (         1:ni, 1:nj, 1:nk))
      allocate(grid%blk(b)%gas_phase%R  (         1:ni, 1:nj, 1:nk))

    end do
    grid%nbound = nbound
  end subroutine Setup_Data_Structure


  subroutine Deallocate_Data(grid)
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    implicit none
    type(GROOT_domain_type), intent(inout) :: grid
    integer :: b

    do b = 1, grid%nb
      deallocate(grid%blk(b)%x,    &
                 grid%blk(b)%y,    &
                 grid%blk(b)%z,    &
                 grid%blk(b)%gas_phase%rho, &
                 grid%blk(b)%gas_phase%T,   &
                 grid%blk(b)%gas_phase%p,   &
                 grid%blk(b)%gas_phase%R)
    end do
  end subroutine Deallocate_Data

end module GROOT_Mod_Allocate_Data
