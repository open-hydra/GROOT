module GROOT_Mod_Bound
  use iso_fortran_env, only: I4 => int32, R8 => real64
  implicit none
  private
  public :: Setup_Bound

contains

  subroutine Setup_Bound(grid)
    use GROOT_Mod_Face, only: get_ind
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use GROOT_Parameters_m, only: MOSE_phase_prefix
    implicit none
    type(GROOT_domain_type), intent(inout), target :: grid
    integer :: nbc
    integer :: i, unitfile, ios
    integer :: j, k, ib, iface, ibc, ind1, ind2
    integer :: ib2, i2, j2, k2, ifa2
    real(R8) :: p, appo

    open(newunit=unitfile,file='INPUT/'//trim(MOSE_phase_prefix)//'bc.txt',status='old',iostat=ios,action='read')
    if (ios /= 0) error stop ' BC file not found'

    nbc = 0
    do
      read(unitfile, *, iostat=ios) ib, i, j, k, iface, ibc
      nbc = nbc + 1
      if (ios /= 0) exit

      call get_ind(i, j, k, iface, ind1, ind2)
      grid%blk(ib)%face(iface)%bc(ind1, ind2) = ibc

      select case (ibc)
        case (101,103)
          read(unitfile, *, iostat=ios) ib2, i2, j2, k2, ifa2
          if (.not. allocated(grid%blk(ib)%face(iface)%bcon)) then
            allocate(grid%blk(ib)%face(iface)%bcon( &
              grid%blk(ib)%face(iface)%n(1), &
              grid%blk(ib)%face(iface)%n(2), 1:5))
          end if
          grid%blk(ib)%face(iface)%bcon(ind1, ind2, 1) = ib2
          grid%blk(ib)%face(iface)%bcon(ind1, ind2, 2) = i2
          grid%blk(ib)%face(iface)%bcon(ind1, ind2, 3) = j2
          grid%blk(ib)%face(iface)%bcon(ind1, ind2, 4) = k2
          grid%blk(ib)%face(iface)%bcon(ind1, ind2, 5) = ifa2
        case(201,401:407, 410, 420, 501:502)
          read(unitfile, *, iostat=ios)
        case(301:309)
          read(unitfile, *, iostat=ios) appo,appo, grid%blk(ib)%face(iface)%eps(ind1,ind2)
        case(102)
          error stop ' BC type 102 not implemented yet'
      end select

      if (ios /= 0) then
        write(*, '(A)') '  ERROR: reading BC file'
        stop
      end if
    end do

    close(unitfile)
  end subroutine Setup_Bound

end module GROOT_Mod_Bound
