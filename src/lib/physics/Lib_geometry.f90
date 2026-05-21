!> @brief This module contains procedures for the manipulation of grid data
module GROOT_Lib_Geometry
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private
  public :: getP0fb, getP0cb, distance, co_area, co_volume, co_normal, checkaxis

contains

!> @brief Suboutine computes cell face center and the rotation matrix
!> from the cell face frame to the universal frame
!>
!> The cell face frame \f$ \mathcal{R}_\mathrm{face}\f$ is identified by its origin
!> \f$ \vec{P}_0 \f$ (face center) and a rotation matrix \f$ R \f$. The rotation 
!> matrix is obtaining expressing the components of the base \f$ (\vec{\tau}, \vec{n}, \vec{b})\f$ in the universal frame
!> For a given cell (i,j,k) the frame is oriented as follows:
!> \image html face_frame.png "Cell face frames definition"
!> and the rotation matrix is obtained as \f$ R = [ \vec{\tau},\vec{n},\vec{b} ] \f$ 
!> While the definition of \f$ \vec{\tau}\f$ and \f$ \vec{b}\f$ is somewhat arbitrary, \f$ \vec{n}\f$ 
!> is chosen to be always the inner facing normal.
subroutine getP0fb(grid, ib, ifa, i1, i2, P0, R, bcm)
  use GROOT_Lib_Math, only: cross, normL2
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  implicit none
  type(GROOT_domain_type),  intent(in) :: grid         !< Block-level data 
  integer,                  intent(in) :: ib,ifa,i1,i2 !< Block, face and cell indexes
  real(kind=R8),            intent(out):: P0(1:3)      !< Face center
  real(kind=R8),            intent(out):: R(1:3,1:3)   !< Rotation matrix
  integer,                  intent(out):: bcm(1:5)     !< Array with cell indexes on block

  integer  :: i, j, k
  ! versors defining face frame
  real(kind=R8) :: tau(1:3), b(1:3), n(1:3)
  real(kind=R8) :: d1(1:3), d2(1:3)

! find block indexes from face indexes
  bcm(1) = ib
  bcm(2) = ifa
  select case (ifa)
    case (1)
      bcm(3) = 1!i
      bcm(4) = i1!j
      bcm(5) = i2!k
    case (2)
      bcm(3) = grid%blk(ib)%dim(1)!i
      bcm(4) = i1!j
      bcm(5) = i2!k
    case (3)
      bcm(3) = i1!i
      bcm(4) = 1!j
      bcm(5) = i2!k
    case (4) 
      bcm(3) = i1!i
      bcm(4) = grid%blk(ib)%dim(2)!j
      bcm(5) = i2!k
    case (5) 
      bcm(3) = i1
      bcm(4) = i2
      bcm(5) = 1  
    case (6)
      bcm(3) = i1
      bcm(4) = i2
      bcm(5) = grid%blk(ib)%dim(3)
  end select
  
  i = bcm(3)
  j = bcm(4)
  k = bcm(5)

! face center coordinates
  P0(1) = 0.25*(grid%blk(ib)%face(ifa)%x(i1,i2)+grid%blk(ib)%face(ifa)%x(i1,i2+1)+ &
                grid%blk(ib)%face(ifa)%x(i1+1,i2)+grid%blk(ib)%face(ifa)%x(i1+1,i2+1))  
  P0(2) = 0.25*(grid%blk(ib)%face(ifa)%y(i1,i2)+grid%blk(ib)%face(ifa)%y(i1,i2+1)+ &
                grid%blk(ib)%face(ifa)%y(i1+1,i2)+grid%blk(ib)%face(ifa)%y(i1+1,i2+1))
  P0(3) = 0.25*(grid%blk(ib)%face(ifa)%z(i1,i2)+grid%blk(ib)%face(ifa)%z(i1,i2+1)+ &
                grid%blk(ib)%face(ifa)%z(i1+1,i2)+grid%blk(ib)%face(ifa)%z(i1+1,i2+1))

! compute face normal as cross product of face diagonals d1 and d2
  d1(1) = grid%blk(ib)%face(ifa)%x(i1+1,i2+1)-grid%blk(ib)%face(ifa)%x(i1,i2)
  d1(2) = grid%blk(ib)%face(ifa)%y(i1+1,i2+1)-grid%blk(ib)%face(ifa)%y(i1,i2)  
  d1(3) = grid%blk(ib)%face(ifa)%z(i1+1,i2+1)-grid%blk(ib)%face(ifa)%z(i1,i2)

  d2(1) = grid%blk(ib)%face(ifa)%x(i1,i2+1) - grid%blk(ib)%face(ifa)%x(i1+1,i2)
  d2(2) = grid%blk(ib)%face(ifa)%y(i1,i2+1) - grid%blk(ib)%face(ifa)%y(i1+1,i2)
  d2(3) = grid%blk(ib)%face(ifa)%z(i1,i2+1) - grid%blk(ib)%face(ifa)%z(i1+1,i2)
! adjuct sign of d2 to ensure that n is the inner facing normal
  select case (ifa)
    case (2)
      d2 = -d2 
    case (3) 
      d2 = -d2
    case (6)
      d2 = -d2
  end select

! compute face normal
  call cross(d1, d2, n)
  n = n/normL2(n)

! compute tau versor directed along the main face direction
  tau(1) = 0.5*(grid%blk(ib)%face(ifa)%x(i1+1,i2)+grid%blk(ib)%face(ifa)%x(i1+1,i2+1)) - &
           0.5*(grid%blk(ib)%face(ifa)%x(i1,i2)+grid%blk(ib)%face(ifa)%x(i1,i2+1))
  tau(2) = 0.5*(grid%blk(ib)%face(ifa)%y(i1+1,i2)+grid%blk(ib)%face(ifa)%y(i1+1,i2+1)) - &
           0.5*(grid%blk(ib)%face(ifa)%y(i1,i2)+grid%blk(ib)%face(ifa)%y(i1,i2+1))
  tau(3) = 0.5*(grid%blk(ib)%face(ifa)%z(i1+1,i2)+grid%blk(ib)%face(ifa)%z(i1+1,i2+1)) - &
           0.5*(grid%blk(ib)%face(ifa)%z(i1,i2)+grid%blk(ib)%face(ifa)%z(i1,i2+1))

  tau = tau/normL2(tau)

! computed b versor to complete the cell face frame
  call cross(n, tau, b)
  b = b/normL2(b)  
! rotation matrix 
  R(1:3,1) = tau(1:3)
  R(1:3,2) = b(1:3)
  R(1:3,3) = n(1:3)

end subroutine getP0fb

!> @brief Subroutine computes cell center and rotation matrix from cell frame to global frame
!> 
!> At the moment the cell frame is taken aligned with the universal frame
subroutine getP0cb(grid, ib, i, j, k, P0, R)
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  implicit none
  type(GROOT_domain_type),  intent(in) :: grid        !< Block-level data
  integer,                  intent(in) :: ib, i, j, k !< Block and cell indexes
  real(kind=R8),            intent(out):: P0(1:3)     !< Cell center
  real(kind=R8),            intent(out):: R(1:3,1:3)  !< Rotation matrix

  ! compute P0 as cell center (average of 8 corner nodes)
  associate(bx => grid%blk(ib)%x, by => grid%blk(ib)%y, bz => grid%blk(ib)%z)
    P0(1) = 0.125*(bx(i,j,k) + bx(i,j+1,k) + bx(i,j,k+1) + bx(i,j+1,k+1) &
                 + bx(i+1,j,k) + bx(i+1,j+1,k) + bx(i+1,j,k+1) + bx(i+1,j+1,k+1))
    P0(2) = 0.125*(by(i,j,k) + by(i,j+1,k) + by(i,j,k+1) + by(i,j+1,k+1) &
                 + by(i+1,j,k) + by(i+1,j+1,k) + by(i+1,j,k+1) + by(i+1,j+1,k+1))
    P0(3) = 0.125*(bz(i,j,k) + bz(i,j+1,k) + bz(i,j,k+1) + bz(i,j+1,k+1) &
                 + bz(i+1,j,k) + bz(i+1,j+1,k) + bz(i+1,j,k+1) + bz(i+1,j+1,k+1))
  end associate
  ! cell frame aligned with universal frame
  R(1:3,1) = (/1.d0, 0.d0, 0.d0/)
  R(1:3,2) = (/0.d0, 1.d0, 0.d0/)
  R(1:3,3) = (/0.d0, 0.d0, 1.d0/)

end subroutine getP0cb

!> @brief Function computes distance between two points 
!>
!> Given point \f$P_0\f$ and \f$P_1\f$ computes \f$ d = \sqrt((P_0(1)-P_1(1))^2+(P_0(2)-P_1(2))^2+(P_0(3)-P_1(3))^2)\f$
function distance(P0, P1) result(dist)
  implicit none
  real(kind=R8), intent(in) :: P0(1:3), P1(1:3) !< Points
  real(kind=R8)             :: dist !< Distance

  dist = sqrt((P0(1) - P1(1))**2 + (P0(2) - P1(2))**2 + (P0(3)- P1(3))**2)
end function distance

!> @brief Function computes cell face area
!>
!> Area is computed as 
!> \f$ A = \frac{1}{2} ||\vec{n}||\f$, where the vector \f$\vec{n}\f$ is computed
!> as the cross product of the face diagonals \f$ \vec{d}_1 \f$ and \f$ \vec{d}_2 \f$
function co_area(face, i1, i2) result(A)
  use GROOT_Mod_Face, only: SBface
  use GROOT_Lib_Math, only: cross, normL2
  use GROOT_Global_m, only: twoDax
  implicit none
  type(SBface),  intent(in) :: face   !< Face-level data
  integer,       intent(in) :: i1, i2 !< Cell indexes on face
  real(kind=R8) :: A
  real(kind=R8) :: d1(1:3),d2(1:3),n(1:3)

  d1(1) = face%x(i1+1,i2+1)-face%x(i1,i2)
  d1(2) = face%y(i1+1,i2+1)-face%y(i1,i2)  
  d1(3) = face%z(i1+1,i2+1)-face%z(i1,i2)

  d2(1) = face%x(i1,i2+1) - face%x(i1+1,i2)
  d2(2) = face%y(i1,i2+1) - face%y(i1+1,i2)
  d2(3) = face%z(i1,i2+1) - face%z(i1+1,i2)

  call cross(d1, d2, n)
  ! area can be computed as 1/2 ||d1 x d2||, which is the same as
  ! ||l1 x l2||, where l1 and l2 are the sides of the quadrilateral and
  ! d the diagonali
  A = 0.5*normL2(n)
  if (twoDax) then
    A = A*360
  endif      

end function co_area


!> @brief Function computes cell volume
function co_volume(blck, i, j, k) result(V)
  use GROOT_Advanced_Types_m, only: GROOT_block_type
  implicit none
  type(GROOT_block_type),  intent(in) :: blck   !< Block-level data
  integer,                 intent(in) :: i, j, k !< Cell indexes in block
  real(kind=R8) :: V
  real(kind=R8), dimension(1:3) :: v1,v2,v3,v4,v5,v6,v7,v8
  real(kind=R8) :: vol1, vol2, vol3, vol4, vol5

  ! cell is split in 5 tetrahedrons
  ! for more info see Dompierre et al, How to subdivide Pyramids, Prisms and
  ! Hexahedra into tetrahedra, 1999

  v1(1) = blck%x(  i,  j,   k); v1(2) = blck%y(  i,   j,   k); v1(3) = blck%z(  i,   j,   k)  
  v2(1) = blck%x(  i,  j, k+1); v2(2) = blck%y(  i,   j, k+1); v2(3) = blck%z(  i,   j, k+1)
  v3(1) = blck%x(i+1,  j, k+1); v3(2) = blck%y(i+1,   j, k+1); v3(3) = blck%z(i+1,   j, k+1)
  v4(1) = blck%x(i+1,  j,   k); v4(2) = blck%y(i+1,   j,   k); v4(3) = blck%z(i+1,   j,   k)
  v5(1) = blck%x(i  ,j+1,   k); v5(2) = blck%y(  i, j+1,   k); v5(3) = blck%z(  i, j+1,   k)
  v6(1) = blck%x(i  ,j+1, k+1); v6(2) = blck%y(  i, j+1, k+1); v6(3) = blck%z(  i, j+1, k+1)
  v7(1) = blck%x(i+1,j+1, k+1); v7(2) = blck%y(i+1, j+1, k+1); v7(3) = blck%z(i+1, j+1, k+1)
  v8(1) = blck%x(i+1,j+1,   k); v8(2) = blck%y(i+1, j+1,   k); v8(3) = blck%z(i+1, j+1,   k)

  vol1 = co_vol_tetra(v1,v2,v3,v6)
  vol2 = co_vol_tetra(v1,v3,v4,v8)
  vol3 = co_vol_tetra(v1,v6,v8,v5)
  vol4 = co_vol_tetra(v3,v6,v8,v7)
  vol5 = co_vol_tetra(v1,v3,v8,v6)

  V = vol1+vol2+vol3+vol4+vol5

end function co_volume

!> @brief Function computes volume of tetrahedron from the vectors specifying vertex positions
function co_vol_tetra(a, b, c, d) result(V)
  use GROOT_Lib_Math, only: cross, scalar
  implicit none
  real(kind=R8), intent(in) :: a(1:3),b(1:3),c(1:3),d(1:3) ! vertexes of tetrahedron
  real(kind=R8) :: V

  real(kind=R8) :: v1(1:3), v2(1:3), v3(1:3), v4(1:3)

  v1 = a-d
  v2 = b-d
  v3 = c-d
  ! V = |(a-d)°((b-d)x(c-d))|/6
  call cross(v2, v3, v4)
  V = abs(scalar(v1,v4))/6
end function co_vol_tetra

!> @brief Subroutine computes the inner facing normal n to face "face" in cell i1,i2
subroutine co_normal(face, ifa, i1, i2, n)
  use GROOT_Mod_Face, only: SBface
  use GROOT_Lib_Math, only: cross, normL2
  implicit none
  type(SBface),  intent(in) :: face   !< Face
  integer,       intent(in) :: ifa    !< Face index
  integer,       intent(in) :: i1, i2 !< Cell indexes on face
  real(kind=R8), intent(out):: n(1:3) !< Normal

  real(kind=R8) :: d1(1:3), d2(1:3)

! compute face normal as cross product of face diagonals d1 and d2
  d1(1) = face%x(i1+1,i2+1)-face%x(i1,i2)
  d1(2) = face%y(i1+1,i2+1)-face%y(i1,i2)  
  d1(3) = face%z(i1+1,i2+1)-face%z(i1,i2)

  d2(1) = face%x(i1,i2+1) - face%x(i1+1,i2)
  d2(2) = face%y(i1,i2+1) - face%y(i1+1,i2)
  d2(3) = face%z(i1,i2+1) - face%z(i1+1,i2)
! adjuct sign of d2 to ensure that n is the inner facing normal
  select case (ifa)
    case (2)
      d2 = -d2 
    case (3) 
      d2 = -d2
    case (6)
      d2 = -d2
  end select
! compute face normal
  call cross(d1, d2, n)
  n = n/normL2(n)
end subroutine co_normal



!> @brief Subroutine check if cell i,j,k is on the symmetry axis
!> warning: only the x axis can be a symmetry axis at the moment
subroutine checkaxis(blck, i, j, k, axis)
  use GROOT_Advanced_Types_m, only: GROOT_block_type
  use GROOT_Global_m, only: twoDax
  implicit none
  type(GROOT_block_type),   intent(in) :: blck
  integer,                  intent(in) :: i, j, k
  logical,                  intent(out):: axis

  axis = .false.
  if ((j==1).and.(blck%face(3)%BC(i,k)==300)) then
    if (twoDax) then
      axis = .true.
    endif
  endif  
  return
end subroutine checkaxis

end module GROOT_Lib_Geometry
