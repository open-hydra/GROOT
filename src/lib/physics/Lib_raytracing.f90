!> @brief This module contains the ray tracing procedures
!> @todo maxiter: allow imax to change if needed
!> @todo index: replace bcm with something nicer
!> @bug  rays: some rays get lost...
module GROOT_Lib_Raytracing
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none
  private 

  public :: getray
  public :: march
contains

!> @brief Subroutine computes ray versor s from azimuth phi, elevation theta and rotation matrix R
subroutine getray(phi, theta, R, s)
  use GROOT_Lib_Math, only: normL2
  implicit none
  real(kind=R8), intent(in) :: phi         !< Azimuth angle
  real(kind=R8), intent(in) :: theta       !< Elevation angle
  real(kind=R8), intent(in) :: R(1:3,1:3)  !< Rotation matrix from cell face frame to universal frame
  real(kind=R8), intent(out):: s(1:3)      !< Ray direction in universal frame

  real(kind=R8)  :: s0(1:3) !< ray direction in cell frame

! compute ray versor in cell face frame  
  s0(1) = cos(theta)*cos(phi)
  s0(2) = cos(theta)*sin(phi)
  s0(3) = sin(theta) 
! rotate versor from cell face frame to universal frame
  s(1) = R(1,1)*s0(1)+R(1,2)*s0(2)+R(1,3)*s0(3)
  s(2) = R(2,1)*s0(1)+R(2,2)*s0(2)+R(2,3)*s0(3)
  s(3) = R(3,1)*s0(1)+R(3,2)*s0(2)+R(3,3)*s0(3)
! normalize (redundant)
 
  s = s/normL2(s)

end subroutine getray


!> @brief Subroutine propagates ray s from P0 until wall is met, then integrates RTE backwads to intensity impinging at the wall
subroutine march(grid, s, P0, bcm, pos, ds, nstep, face, rayflag)
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  use GROOT_Lib_Geometry, only: distance
  use GROOT_Global_m, only: imax
  implicit none
  type(GROOT_domain_type),  intent(in)    :: grid             !< Block-level data 
  real(kind=R8),            intent(inout) :: s(1:3)           !< Ray versor            
  real(kind=R8),            intent(in)    :: P0(1:3)          !< Starting point
  integer,                  intent(in)    :: bcm(1:5)         !< Array containing the cell and face indexes of P0
  integer,                  intent(out)   :: pos(1:imax,1:4)  !< Array with cell and block index
  real(kind=R8),            intent(out)   :: ds(1:imax)       !< Array width distances travelled by ray
  integer,                  intent(out)   :: nstep            !< Number of steps
  integer,                  intent(out)   :: face             !< Face id (1-6) 
  logical,                  intent(out)   :: rayflag          !< Error flag

  real(kind=R8)  :: Pend(1:3),Pstart(1:3) !< start/end points
  integer   :: i, j, k, ib, ii, ifa, ifa1, ifa2, i2, j2, k2, ib2
  logical   :: vai
  
  rayflag = .true.
! (i, j, k, ifa, ib) of starting point are known
  ib  = bcm(1)
  ifa = bcm(2)
  i   = bcm(3)
  j   = bcm(4) 
  k   = bcm(5) 
 
  vai = .true.
  ii = 1
  
  Pstart = P0

  do while (vai)
    ! find intersection Pend of line P0+s*t and cell walls
    call propagate(grid, Pstart, i, j, k, ib, ifa, s, Pend, ifa1, rayflag)
    ! compute ds, je, k
    if (.not.rayflag) then
      exit
    endif
    ds(ii)    = distance(Pstart,Pend)
    pos(ii,1) = ib; pos(ii,2)=i; pos(ii,3)=j; pos(ii,4)=k
    ! find (i,j,k,ib) of next cell  
    call getnext(grid, i, j, k, ib, ifa1, i2, j2, k2, ib2, ifa2, vai, s, Pend)

    ib  = ib2
    i   = i2
    j   = j2
    k   = k2
    ifa = ifa2
    Pstart = Pend
    ii = ii + 1
    if (ii>=imax) then
      print *, 'too many steps?'
      stop
    endif
  enddo

  face = ifa1
  nstep = ii-1
end subroutine march

!> @brief Subroutine finds intersection between ray s starting from P0 and the cell faces
subroutine propagate(grid, P0, i, j, k, ib, ifa, s, P, ifa1, rayflag)
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  use GROOT_Lib_Geometry, only: checkaxis
  implicit none
  type(GROOT_domain_type),  intent(in)    :: grid              !< Block-level data
  real(kind=R8),            intent(in)    :: P0(1:3)           !< Starting point
  real(kind=R8),            intent(in)    :: s(1:3)            !< Ray versor
  integer,                  intent(in)    :: i, j, k, ib, ifa  !< Cell indexes and face where P0 is
  real(kind=R8),            intent(out)   :: P(1:3)            !< Point of intersection btw s and cell boundary
  integer,                  intent(out)   :: ifa1              !< Face where P is
  logical,                  intent(out)   :: rayflag           !< Error flag

  real(kind=R8)  :: M(1:3,1:4)
  integer   :: iface, geom
  logical   :: flag, axis 

  rayflag = .true.

  ! in 3D computations cell faces are planes, specified through the M matrix
  do iface = 1, 6
    ! M is a 3x4 matrix with the vertex coordinates
    ! vertex order is abritrarily chosen anti-clockwise

    ! check if cell i,j,k is on the symmetry axis
    call checkaxis(blck = grid%blk(ib), i = i, j = j, k = k, axis =  axis)
    geom = 4
    if (iface.ne.ifa) then
      associate(bx => grid%blk(ib)%x, by => grid%blk(ib)%y, bz => grid%blk(ib)%z)
      select case (iface)
        case (1)
          M(1,:) = [bx(i,j,k+1), bx(i,j,k), bx(i,j+1,k), bx(i,j+1,k+1)]
          M(2,:) = [by(i,j,k+1), by(i,j,k), by(i,j+1,k), by(i,j+1,k+1)]
          M(3,:) = [bz(i,j,k+1), bz(i,j,k), bz(i,j+1,k), bz(i,j+1,k+1)]
          if (axis) then
            geom = 3
          endif
        case (2)
          M(1,:) = [bx(i+1,j,k+1), bx(i+1,j,k), bx(i+1,j+1,k), bx(i+1,j+1,k+1)]
          M(2,:) = [by(i+1,j,k+1), by(i+1,j,k), by(i+1,j+1,k), by(i+1,j+1,k+1)]
          M(3,:) = [bz(i+1,j,k+1), bz(i+1,j,k), bz(i+1,j+1,k), bz(i+1,j+1,k+1)]
          if (axis) then
            geom = 3
          endif
        case (3)
          M(1,:) = [bx(i+1,j,k), bx(i,j,k), bx(i,j,k+1), bx(i+1,j,k+1)]
          M(2,:) = [by(i+1,j,k), by(i,j,k), by(i,j,k+1), by(i+1,j,k+1)]
          M(3,:) = [bz(i+1,j,k), bz(i,j,k), bz(i,j,k+1), bz(i+1,j,k+1)]
        case (4)
          M(1,:) = [bx(i+1,j+1,k), bx(i,j+1,k), bx(i,j+1,k+1), bx(i+1,j+1,k+1)]
          M(2,:) = [by(i+1,j+1,k), by(i,j+1,k), by(i,j+1,k+1), by(i+1,j+1,k+1)]
          M(3,:) = [bz(i+1,j+1,k), bz(i,j+1,k), bz(i,j+1,k+1), bz(i+1,j+1,k+1)]
        case (5)
          M(1,:) = [bx(i,j,k), bx(i+1,j,k), bx(i+1,j+1,k), bx(i,j+1,k)]
          M(2,:) = [by(i,j,k), by(i+1,j,k), by(i+1,j+1,k), by(i,j+1,k)]
          M(3,:) = [bz(i,j,k), bz(i+1,j,k), bz(i+1,j+1,k), bz(i,j+1,k)]
        case (6)
          M(1,:) = [bx(i,j,k+1), bx(i+1,j,k+1), bx(i+1,j+1,k+1), bx(i,j+1,k+1)]
          M(2,:) = [by(i,j,k+1), by(i+1,j,k+1), by(i+1,j+1,k+1), by(i,j+1,k+1)]
          M(3,:) = [bz(i,j,k+1), bz(i+1,j,k+1), bz(i+1,j+1,k+1), bz(i,j+1,k+1)]
      end select
      end associate

      if ((iface.eq.3).and.(axis)) then
        flag = .false. 
      else  
        call intface3D(P0, s, M, geom, P, flag)   
      endif

      if (flag) then
        ifa1 = iface
        return
      endif

    endif
  enddo
  
  if (axis) then
    call intaxis(blck=grid%blk(ib), i=i,j=j,k=k,P0=P0,s=s,P=P,flag=flag)
    if (flag) then
      ifa1 = 3
      return
    else

    endif
  endif
  rayflag = .false.


end subroutine propagate

!> @brief Subroutine finds intersection between P0+s*t and plane given by the points stored in M
subroutine intface3D(P0, s, M, geom, P, flag)
  use GROOT_Lib_Math, only: normL2, cross, scalar
  use GROOT_Global_m, only: pi
  implicit none
  real(kind=R8),  intent(in) :: P0(1:3)     !< Ray starting point 
  real(kind=R8),  intent(in) :: s(1:3)      !< Ray versor
  real(kind=R8),  intent(in) :: M(1:3,1:4)  !< 3x4 matrix containg the coordinates of the face vertexes
  integer,        intent(in) :: geom        !< identifies axis
  real(kind=R8),  intent(out):: P(1:3)      !< Intersection point
  logical,        intent(out):: flag        !< Error flag

  real(kind=R8) :: d1(1:3),d2(1:3),n(1:3),Pcc(1:3)
  real(kind=R8) :: v1(1:3),v2(1:3), P2(1:3)
  real(kind=R8) :: t,ns,a12,a23,a34,a41,atot,tol

  tol = 1e-7

  flag = .false.
  P = (/0.d0, 0.d0, 0.d0/)
  
  ! normal to the face plane
  d1 = M(1:3,3)-M(1:3,1)
  d2 = M(1:3,4)-M(1:3,2)
  call cross(d1, d2, n)
  n = n/normL2(n)

  ! face center  
  Pcc(1) = 0.25*(M(1,1)+M(1,2)+M(1,3)+M(1,4))
  Pcc(2) = 0.25*(M(2,1)+M(2,2)+M(2,3)+M(2,4))
  Pcc(3) = 0.25*(M(3,1)+M(3,2)+M(3,3)+M(3,4))

  ! intersection between face plane and ray
  ! plane equation is given by n x (P-Pcc) = 0
  ! where n is the plane normal, P is [x,y,z] and Pcc a point on the plane
  ! substituting P = P0+s*t one gets the formula below
  ns = scalar(n, s)
  if (ns.ne.0.d0) then
    P2 = Pcc-P0
    t = scalar(n, P2)/ns
    P = P0 + s*t
    if (t>=0) then
      ! check if sum of angles btw vectors from P to vertexes
      ! is 2*pi, if less P is not inside the polygon
      ! taken from eecs.umich.edu/courses/eecs380/HANDOUTS/PROJ2/InsidePoly.html
      if (geom==4) then
      ! vertexes 1-2
        v1 = M(1:3,1)-P
        v2 = M(1:3,2)-P
        a12 = acos(scalar(v1,v2)/(normL2(v1)*normL2(v2)))
        if (a12/=a12) a12 = 0.d0
      elseif (geom==3) then
        a12 = 0.d0
      endif
      ! vertexes 2-3
      v1 = M(1:3,2)-P
      v2 = M(1:3,3)-P
      a23 = acos(scalar(v1,v2)/(normL2(v1)*normL2(v2)))
      if (a23/=a23) a23 = 0.d0
     ! vertexes 3-4
      v1 = M(1:3,3)-P
      v2 = M(1:3,4)-P
      a34 = acos(scalar(v1,v2)/(normL2(v1)*normL2(v2)))
      if (a34/=a34) a34 = 0.d0
      ! vertexes 4-1
      v1 = M(1:3,4)-P
      v2 = M(1:3,1)-P
      if (geom==3) v2 = M(1:3,2)-P
      a41 = acos(scalar(v1,v2)/(normL2(v1)*normL2(v2)))
      if (a41/=a41) a41 = 0.d0

      atot = a12+a23+a34+a41
     if (atot.ge.(2*pi-tol)) then
        flag = .true.
        return
      endif 
      
    endif
  endif

end subroutine intface3D

!> @brief Subrouitne looks for intersection with the symmetry axis
subroutine intaxis(blck, i, j, k, P0, s, P, flag)
  use GROOT_Advanced_Types_m, only: GROOT_block_type
  use GROOT_Lib_Math, only: scalar
  implicit none
  type(GROOT_block_type),   intent(in) :: blck
  integer,                  intent(in) :: i,j,k
  real(kind=R8),            intent(in) :: P0(1:3)
  real(kind=R8),            intent(in) :: s(1:3)
  real(kind=R8),            intent(out):: P(1:3)
  logical,                  intent(out):: flag

  real(kind=R8) :: n(1:3), t, ns, Pcc(1:3), P2(1:3)

  n = (/0.d0, 1.d0, 0.d0/)
  P = (/-1.d0,-1.d0,-1.d0/)
  Pcc(1) = 0.5*(blck%x(i,j,k)+blck%x(i+1,j,k))
  Pcc(2) = 0.5*(blck%y(i,j,k)+blck%y(i+1,j,k))
  Pcc(3) = 0.d0
  ns = scalar(n, s)
  if (ns.ne.0d0) then
    P2 = Pcc-P0
    t = scalar(n, P2)/ns
    P = P0 + s*t
    if (t.gt.0) then
      flag = .true.
      return
    endif
  endif
end subroutine intaxis

!> @brief Subroutine finds indexes of next cell traversed by the ray
!> if the block boundaries are reached, either a new block is encountered
!> (if BC==1) or the integration is stopped
subroutine getnext(grid, i, j, k, ib, ifa1, i2, j2, k2, ib2, ifa2, vai, s, P)
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  implicit none
  type(GROOT_domain_type),  intent(in)    :: grid                  !< Block-level data
  integer,                  intent(in)    :: i, j, k, ib, ifa1     !< indexes of current cell
  integer,                  intent(out)   :: i2, j2, k2, ib2, ifa2 !< indexes of next cell
  logical,                  intent(out)   :: vai                   !< flag to continue/stop integration
  real(kind=R8),            intent(inout) :: s(1:3)                !< Ray versor
  real(kind=R8),            intent(inout) :: P(1:3)                !< Impinging point
! default values
  vai = .true.
  ib2 = ib
  i2 = i
  j2 = j
  k2 = k
  ifa2 = ifa1
  s = s

  select case (ifa1)

    case (1)
    if (i.ne.1) then
      ib2  = ib 
      i2   = i-1
      j2   = j
      k2   = k
      ifa2 = 2
    else
      select case (grid%blk(ib)%face(1)%bc(j,k))
        case (101,103)
        ! connection b.c
        ib2  = grid%blk(ib)%face(1)%Bcon(j,k,1)
        i2   = grid%blk(ib)%face(1)%bcon(j,k,2)
        j2   = grid%blk(ib)%face(1)%bcon(j,k,3)
        k2   = grid%blk(ib)%face(1)%bcon(j,k,4)
        ifa2 = grid%blk(ib)%face(1)%Bcon(j,k,5)
        case (200)
        ! periodic b.c
        case (300)
          ! symmetry b.b
          call specular(s = s, face = grid%blk(ib)%face(1), ifa = 1, i1 = j, i2 = k)
        case (301:)
          vai = .false.
      end select
    endif

    case (2)
    if (i.ne.grid%blk(ib)%dim(1)) then
      ib2  = ib
      i2   = i+1
      j2   = j
      k2   = k
      ifa2 = 1
    else
      select case (grid%blk(ib)%face(2)%bc(j,k))
        case (101,103)
         ! connection b.c
        ib2  = grid%blk(ib)%face(2)%Bcon(j,k,1)
        i2   = grid%blk(ib)%face(2)%bcon(j,k,2)
        j2   = grid%blk(ib)%face(2)%bcon(j,k,3)
        k2   = grid%blk(ib)%face(2)%bcon(j,k,4)
        ifa2 = grid%blk(ib)%face(2)%Bcon(j,k,5)
        case (200)
          ! periodic b.c

        case (300)
          ! symmetry b.c
          call specular(s = s, face = grid%blk(ib)%face(2), ifa = 2, i1 = j, i2 = k)
        case (301:)
          vai = .false.
      end select
    endif

    case (3)
    if (j.ne.1) then
      ib2  = ib
      i2   = i
      j2   = j-1
      k2   = k
      ifa2 = 4
    else
      select case (grid%blk(ib)%face(3)%bc(i,k))
        case (101,103)
        ib2  = grid%blk(ib)%face(3)%Bcon(i,k,1)
        i2   = grid%blk(ib)%face(3)%bcon(i,k,2)
        j2   = grid%blk(ib)%face(3)%bcon(i,k,3)
        k2   = grid%blk(ib)%face(3)%bcon(i,k,4)
        ifa2 = grid%blk(ib)%face(3)%Bcon(i,k,5)
        case (200)
      
        case (300)
          call specular(s = s, face = grid%blk(ib)%face(3), ifa = 3, i1 = i, i2 = k)
        case (301:)
          vai = .false.
      end select
    endif

    case (4)
    if (j.ne.grid%blk(ib)%dim(2)) then
      ib2  = ib
      i2   = i
      j2   = j+1
      k2   = k
      ifa2 = 3
    else
      select case (grid%blk(ib)%face(4)%bc(i,k))
        case (101,103)
        ib2  = grid%blk(ib)%face(4)%Bcon(i,k,1)
        i2   = grid%blk(ib)%face(4)%bcon(i,k,2)
        j2   = grid%blk(ib)%face(4)%bcon(i,k,3)
        k2   = grid%blk(ib)%face(4)%bcon(i,k,4)
        ifa2 = grid%blk(ib)%face(4)%Bcon(i,k,5)
        case (200)
          ! periodic b.c
        case (300)
          call specular(s = s, face = grid%blk(ib)%face(4), ifa = 4, i1 = i, i2 = k)
        case (301:)
          vai = .false.
      end select
    endif

    case (5)
    if (k.ne.1) then
      ib2  = ib
      i2   = i
      j2   = j
      k2   = k-1
      ifa2 = 6
    else
      select case (grid%blk(ib)%face(5)%bc(i,j))
        case (101,103)
        ib2  = grid%blk(ib)%face(5)%Bcon(i,j,1)
        i2   = grid%blk(ib)%face(5)%Bcon(i,j,2)
        j2   = grid%blk(ib)%face(5)%Bcon(i,j,3)
        k2   = grid%blk(ib)%face(5)%Bcon(i,j,4)
        ifa2 = grid%blk(ib)%face(5)%Bcon(i,j,5)
        case (200)
          call periodic(grid, ib, 5, i, j, P, s, ib2, i2, j2, k2, ifa2)
        case (300)
          call specular(s = s, face = grid%blk(ib)%face(5), ifa = 5, i1 = i, i2 = j)
        case (301:)
          vai = .false.
      end select
    endif
 
    case (6)
    if (k.ne.grid%blk(ib)%dim(3)) then
      ib2  = ib
      i2   = i
      j2   = j
      k2   = k+1
      ifa2 = 5
    else
      select case (grid%blk(ib)%face(6)%bc(i,j))
        case (101,103)
        ib2  = grid%blk(ib)%face(6)%Bcon(i,j,1)
        i2   = grid%blk(ib)%face(6)%Bcon(i,j,2)
        j2   = grid%blk(ib)%face(6)%Bcon(i,j,3)
        k2   = grid%blk(ib)%face(6)%Bcon(i,j,4)
        ifa2 = grid%blk(ib)%face(6)%Bcon(i,j,5)
        case (200)
          call periodic(grid, ib, 6, i, j, P, s, ib2, i2, j2, k2, ifa2)
        case (300)
          call specular(s = s, face = grid%blk(ib)%face(6), ifa = 6, i1 = i, i2 = j)
        case (301:)
          vai = .false.
      end select
    endif 

  end select
 
end subroutine getnext

!> @brief Function computes the scalar product of (P1-P2) and s
function cosine(P1,P2,s) result(c)
  implicit none
  real(kind=R8), intent(in) :: P1(1:3),P2(1:3)
  real(kind=R8), intent(in) :: s(1:3)
  real(kind=R8)             :: c

  c = (P1(1)-P2(1))*s(1) + (P1(2)-P2(2))*s(2) + (P1(3)-P2(3))*s(3)

end function cosine

!> @brief Subroutine reflects the ray s with respect to the symmetry surface face
subroutine specular(s, face, ifa, i1, i2)
  use GROOT_Mod_Face, only: SBface 
  use GROOT_Lib_Geometry, only: co_normal 
  use GROOT_Lib_Math, only: scalar
  implicit none
  real(kind=R8),      intent(inout) :: s(1:3)  !< Ray versor
  type(SBface),  intent(in)    :: face    !< Symmetry face
  integer,       intent(in)    :: ifa     !< Face index
  integer,       intent(in)    :: i1, i2  !< Indexes on symmetry face

  real(kind=R8) :: n(1:3)
  ! compute the face normal
  call co_normal(face, ifa, i1, i2, n)
  ! reflect the ray s
  ! s = s_perp + s_par, s_par = (s°n)n, s_per = s-s_par
  ! s_ref = s_perp - s_par, change the sign of the normal component
  s = s-2*scalar(s,n)*n

end subroutine specular

!> @brief Subroutine applies a periodic boundary condition to the ray s 
!> @note this function works only for faces 5 and 6 of the same block
subroutine periodic(grid, ib_1, ifa_1, i1_1, i2_1, P_1, s, ib_2, i_2, j_2, k_2, ifa_2)
  use GROOT_Lib_Math, only: mat_mul_vec, invert 
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  use GROOT_Lib_Geometry, only: getP0fb
  implicit none
  type(GROOT_domain_type),  intent(in)    :: grid                   !< Block-level data
  integer,                  intent(in)    :: ib_1, ifa_1, i1_1, i2_1    !< Indexes of "inlet" face
  real(kind=R8),            intent(inout) :: P_1(1:3)                   !< Impinging point on "inlet" face
  real(kind=R8),            intent(inout) :: s(1:3)                     !< Ray versor
  integer,                  intent(out)   :: ib_2, ifa_2, i_2, j_2, k_2 !< Indexes of "outlet" face

  real(kind=R8) :: R_1(1:3,1:3), R_2(1:3,1:3)
  real(kind=R8) :: P0_1(1:3), P0_2(1:3)
  integer       :: bcm_1(1:5), bcm_2(1:5), i1_2, i2_2
  real(kind=R8) :: s_loc1(1:3), s_loc2(1:3)
  real(kind=R8) :: P_loc1(1:3), P_loc2(1:3)
  real(kind=R8) :: invR_1(1:3,1:3), P1diff(1:3)
  ib_2 = ib_1    

  ! it is assumed that the two connected faces have the same discretization
  ! namely i,j are the same 
  select case (ifa_1) 
    case(1,2,3,4)
      print *, 'periodic boundary condition not implemented'
      stop
    case (5)
      ifa_2 = 6      
      k_2 = grid%blk(ib_2)%dim(3)
      j_2 = i2_1
      i_2 = i1_1
    case (6)
      ifa_2 = 5
      k_2 = 1
      j_2 = i2_1
      i_2 = i1_1
  end select
  i1_2 = i1_1
  i2_2 = i2_1
  ! get center and rotation matrix of face frames
  call getP0fb(grid, ib_1, ifa_1, i1_1, i2_1, P0_1, R_1, bcm_1)
  ! beware: R yields a local frame which has always inner facing normal
  call getP0fb(grid, ib_2, ifa_2, i1_2, i2_2, P0_2, R_2, bcm_2)

  ! rotate ray
  ! first we get the ray in the local reference frame of the inlet face, s_loc1 = inv(R_1)*s
  ! then we note that, in the local frame of the outlet face, the ray has opposed direction s_loc2 = -s_loc1
  ! this is because the frames are defined in order to ensure that the normal n is always inward
  ! finally we compute the components of the rotated ray in the global frame s = R_2*s_loc2 
  invR_1 = invert(R_1)  
  s_loc1 = mat_mul_vec(invR_1,s)
  s_loc2 = s_loc1 ; s_loc2(3) = -s_loc2(3) ; s_loc2(2) = -s_loc2(2)
  s = mat_mul_vec(R_2,s_loc2)

  ! traslate application point
  ! first we compute the position of the impinging point in the local frame of the inlet face, P_loc1 = inv(R1)*(P1-P01)
  ! due to the different orientation of the two local reference frames, to obtain the coordinates of the outlet point
  ! on the outlet face we have to reverse the sign of second component P_loc2 = P_loc1*[1,-1,1]
  ! we can than compute the coordinates of the outlet point in the global frame

  P1diff = P_1-P0_1
  P_loc1 = mat_mul_vec(invR_1,P1diff)
  P_loc2 = P_loc1 ; P_loc2(2) = -P_loc2(2)
  P_1 = P0_2 + mat_mul_vec(R_2,P_loc2)
end subroutine periodic


end module GROOT_Lib_Raytracing

