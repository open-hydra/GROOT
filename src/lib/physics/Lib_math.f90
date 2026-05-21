!> @brief This module contains mathematical procedures for the manipulation of vectors and tensors

module GROOT_Lib_Math
  use, intrinsic :: iso_fortran_env, only : I4 => int32, R8 => real64
  implicit none

  contains

  !> @brief Function for the computation of the determinant of a 3x3 matrix 
  function determinant(ten) result(det)
  implicit none
  real(kind=R8), intent(in):: ten(1:3,1:3) !< Tensor.
  real(kind=R8)            :: det          !< Determinant.
  det= ten(1,1)*(ten(3,3)*ten(2,2)-ten(3,2)*ten(2,3))&
      -ten(2,1)*(ten(3,3)*ten(1,2)-ten(3,2)*ten(1,3))&
      +ten(3,1)*(ten(2,3)*ten(1,2)-ten(2,2)*ten(1,3))
  return
  end function determinant

  !> Function for 3x3 matrix inversion 
  function invert(ten) result(inv)
  implicit none
  real(kind=R8), intent(in)  :: ten(1:3,1:3) !< Tensor to be inverted.
  real(kind=R8)              :: inv(1:3,1:3) !< Tensor inverted.
  real(kind=R8)              :: det          !< Determinant and 1/Determinant.
  det = determinant(ten)
  if (det/=0.d0) then
    det = 1.d0/det
    inv(1,1) = (ten(2,2)*ten(3,3)-ten(2,3)*ten(3,2))*det
    inv(2,1) = -(ten(2,1)*ten(3,3)-ten(2,3)*ten(3,1))*det
    inv(3,1) = (ten(2,1)*ten(3,2)-ten(2,2)*ten(3,1))*det
    inv(1,2) = -(ten(1,2)*ten(3,3)-ten(1,3)*ten(3,2))*det
    inv(2,2) = (ten(1,1)*ten(3,3)-ten(1,3)*ten(3,1))*det
    inv(3,2) = -(ten(1,1)*ten(3,2)-ten(1,2)*ten(3,1))*det
    inv(1,3) = (ten(1,2)*ten(2,3)-ten(1,3)*ten(2,2))*det
    inv(2,3) = -(ten(1,1)*ten(2,3)-ten(1,3)*ten(2,1))*det
    inv(3,3) = (ten(1,1)*ten(2,2)-ten(1,2)*ten(2,1))*det
  else
    print *, 'matrix is not invertible'
    stop     
  endif
  
  return
  endfunction invert

  !> @brief Function for matrix-vector multiplication
  function mat_mul_vec(ten, vec) result(vec2)
    implicit none
    real(kind=R8), intent(in) :: ten(1:3,1:3)
    real(kind=R8), intent(in) :: vec(1:3)
    real(kind=R8) :: vec2(1:3)

    vec2(1) = ten(1,1)*vec(1)+ten(1,2)*vec(2)+ten(1,3)*vec(3)
    vec2(2) = ten(2,1)*vec(1)+ten(2,2)*vec(2)+ten(2,3)*vec(3)
    vec2(3) = ten(3,1)*vec(1)+ten(3,2)*vec(2)+ten(3,3)*vec(3)

  end function mat_mul_vec


!> @brief Function computes the L2 norm of a 3-component vector
function normL2(v) result(norm)
  implicit none
  real(kind=R8), intent(in) :: v(1:3) !< Vector
  real(kind=R8) :: norm

  norm = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
end function normL2

!> @brief Function computes cross product between two 3-component vectors
subroutine cross(v1, v2, v)
  implicit none
  real(kind=R8), intent(in) :: v1(1:3)
  real(kind=R8), intent(in) :: v2(1:3)
  real(kind=R8), intent(out):: v(1:3)

  v(1) = (v1(2)*v2(3)-v1(3)*v2(2))
  v(2) = (v1(3)*v2(1)-v1(1)*v2(3))
  v(3) = (v1(1)*v2(2)-v1(2)*v2(1))

end subroutine cross

!> @brief Function computes scalar product between two 3-component vectors 
function scalar(v1, v2) result(s)
  implicit none
  real(kind=R8), intent(in) :: v1(1:3), v2(1:3)
  real(kind=R8) :: s
  
  s = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)

end function scalar

end module GROOT_Lib_Math
