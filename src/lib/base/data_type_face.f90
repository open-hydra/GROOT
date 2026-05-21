!! @brief Definition and procedures for the SBface type (structured-block face).
module GROOT_Mod_Face
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use GROOT_Global_m, only: Ngg

  implicit none
  private

  public :: SBface
  public :: get_ind       !< Face indexes from block indexes
  public :: set_bc        !< Assign a single BC value to all face cells
  public :: face_to_block !< Block indexes from face-surface indexes

  !! ------------------------------------------------------
  type :: SBface
    integer :: n(2)                                     !< Cell count on face
    real(R8), allocatable :: x(:,:), y(:,:), z(:,:)    !< Nodal coordinates
    real(R8), allocatable :: T(:,:)                     !< Wall temperature
    real(R8), allocatable :: Q(:,:,:)                   !< Spectral radiative heat flux
    real(R8), allocatable :: Qtot(:,:)                  !< Total radiative heat flux
    real(R8), allocatable :: eps(:,:)                   !< Emissivity
    real(R8), allocatable :: a(:,:,:)                   !< WSGG gray-gas weights
    integer,  allocatable :: bc(:,:)                    !< BC flag per cell
    integer,  allocatable :: bcon(:,:,:)                !< Connected block info (BC type 1)
  contains
    procedure :: alloc => alloc_face
    procedure :: free  => free_face
  end type SBface
  !! ------------------------------------------------------

contains

  subroutine alloc_face(face)
    implicit none
    class(SBface), intent(inout) :: face
    integer :: n1, n2

    n1 = face%n(1)
    n2 = face%n(2)
    call free_face(face)
    allocate( face%x(1:n1+1, 1:n2+1) )
    allocate( face%y(1:n1+1, 1:n2+1) )
    allocate( face%z(1:n1+1, 1:n2+1) )
    allocate( face%T(1:n1, 1:n2) )
    allocate( face%Q(1:n1, 1:n2, 0:Ngg) )
    allocate( face%Qtot(1:n1, 1:n2) )
    allocate( face%a(1:n1, 1:n2, 0:Ngg) )
    allocate( face%eps(1:n1, 1:n2) )
    allocate( face%bc(1:n1, 1:n2) )
  end subroutine alloc_face


  subroutine free_face(face)
    implicit none
    class(SBface), intent(inout) :: face

    if (allocated(face%x))    deallocate(face%x)
    if (allocated(face%y))    deallocate(face%y)
    if (allocated(face%z))    deallocate(face%z)
    if (allocated(face%T))    deallocate(face%T)
    if (allocated(face%Q))    deallocate(face%Q)
    if (allocated(face%Qtot)) deallocate(face%Qtot)
    if (allocated(face%a))    deallocate(face%a)
    if (allocated(face%eps))  deallocate(face%eps)
    if (allocated(face%bc))   deallocate(face%bc)
  end subroutine free_face


  !! @brief Return face-local (i1,i2) from block (ix,iy,iz) and face index.
  subroutine get_ind(ix, iy, iz, iface, i1, i2)
    implicit none
    integer, intent(in)  :: ix, iy, iz, iface
    integer, intent(out) :: i1, i2

    select case (iface)
      case (1, 2);  i1 = iy ; i2 = iz
      case (3, 4);  i1 = ix ; i2 = iz
      case (5, 6);  i1 = ix ; i2 = iy
    end select
  end subroutine get_ind


  !! @brief Return block (i,j,k) and ghost (ii,jj,kk) from face (ind1,ind2).
  subroutine face_to_block(ind1, ind2, iface, nx, ny, nz, i, j, k, ii, jj, kk)
    implicit none
    integer, intent(in)  :: ind1, ind2, iface, nx, ny, nz
    integer, intent(out) :: i, j, k, ii, jj, kk

    select case (iface)
      case (1);  i=1;   j=ind1; k=ind2; ii=0;    jj=ind1; kk=ind2
      case (2);  i=nx;  j=ind1; k=ind2; ii=nx+1; jj=ind1; kk=ind2
      case (3);  i=ind1; j=1;   k=ind2; ii=ind1; jj=0;    kk=ind2
      case (4);  i=ind1; j=ny;  k=ind2; ii=ind1; jj=ny+1; kk=ind2
      case (5);  i=ind1; j=ind2; k=1;   ii=ind1; jj=ind2; kk=0
      case (6);  i=ind1; j=ind2; k=nz;  ii=ind1; jj=ind2; kk=nz+1
    end select
  end subroutine face_to_block


  !! @brief Assign the same BC value to every cell on a face.
  subroutine set_bc(face, bcval)
    implicit none
    class(SBface), intent(inout) :: face
    integer,       intent(in)    :: bcval
    integer :: i1, i2

    do i2 = 1, face%n(2)
      do i1 = 1, face%n(1)
        face%bc(i1, i2) = bcval
      end do
    end do
  end subroutine set_bc

end module GROOT_Mod_Face
