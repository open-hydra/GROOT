module GROOT_Backend_INI
  implicit none
  private
  public :: Load_Ini

contains

  !! Read values from an already-loaded FiNeR ini file into all
  !! registry-pointed variables.  Must be called after all Register_*
  !! routines have populated the registry.
  subroutine Load_Ini(fini)
    use GROOT_Input_Registry, only: reg
    use finer,                only: file_ini
    implicit none
    type(file_ini), intent(in) :: fini
    integer :: i, error

    do i = 1, reg%size
      error = 1

      if (associated(reg%params(i)%value%i)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, &
                      val=reg%params(i)%value%i, error=error)

      else if (associated(reg%params(i)%value%r)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, &
                      val=reg%params(i)%value%r, error=error)

      else if (associated(reg%params(i)%value%l)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, &
                      val=reg%params(i)%value%l, error=error)

      else if (associated(reg%params(i)%value%s)) then
        call fini%get(reg%params(i)%section, reg%params(i)%name, &
                      val=reg%params(i)%value%s, error=error)
      end if

      if (error == 0) reg%params(i)%is_set = .true.
    end do

  end subroutine Load_Ini

end module GROOT_Backend_INI
