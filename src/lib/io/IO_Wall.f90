module GROOT_IO_Wall
  use GROOT_Advanced_Types_m, only: GROOT_domain_type
  use Lib_ORION_Data
  use GROOT_Global_m

  implicit none
  private
  public :: Initialize_Wall_File
  public :: Write_Wall_Solution
  public :: viscous_flag
  public :: IOwall
  public :: Owallnames

  logical, allocatable :: viscous_flag(:,:)   !< (block, face) viscous BC flag

  integer          :: Onwall
  type(ORION_data) :: IOwall
  character(llen)  :: Owallnames

contains

  subroutine Initialize_Wall_File(grid)
    use IR_Precision, only: str
    implicit none
    type(GROOT_domain_type), intent(in) :: grid
    integer :: b, f, bb, i, j, k

    Onwall    = 2
    Owallnames = '"emissivity" "Q"'
    if (trim(model) == 'wsgg-H2O' .or. trim(model) == 'wsgg-H2OCO2') then
      Owallnames = trim(Owallnames)//' "a(0)","a(1)","a(2)","a(3)","a(4)"'
      Onwall     = Onwall + 5
    end if

    if (allocated(IOwall%block)) deallocate(IOwall%block)
    allocate(IOwall%block(1:count(viscous_flag)))

    bb = 0
    do b = 1, grid%nb
      do f = 1, 6
        if (.not. viscous_flag(b, f)) cycle
        bb = bb + 1
        IOwall%block(bb)%name = 'B'//trim(str(.true.,b))//'F'//trim(str(.true.,f))

        associate(Nx => grid%blk(b)%dim(1), &
                  Ny => grid%blk(b)%dim(2), &
                  Nz => grid%blk(b)%dim(3))
          select case (f)
            case (1, 2)
              allocate(IOwall%block(bb)%mesh(1:3, 0:0,  0:Ny, 0:Nz))
              allocate(IOwall%block(bb)%vars(1:Onwall, 1, 1:Ny, 1:Nz))
              do k = 0, Nz; do j = 0, Ny
                IOwall%block(bb)%mesh(1, 0, j, k) = grid%blk(b)%x(merge(1,Nx+1,f==1), j+1, k+1)
                IOwall%block(bb)%mesh(2, 0, j, k) = grid%blk(b)%y(merge(1,Nx+1,f==1), j+1, k+1)
                IOwall%block(bb)%mesh(3, 0, j, k) = grid%blk(b)%z(merge(1,Nx+1,f==1), j+1, k+1)
              end do; end do
              IOwall%block(bb)%vars(1, 1, 1:Ny, 1:Nz) = &
                grid%blk(b)%face(f)%eps(1:Ny, 1:Nz)
              IOwall%block(bb)%vars(2, 1, 1:Ny, 1:Nz) = &
                grid%blk(b)%face(f)%Qtot * grid%blk(b)%face(f)%eps &
                - grid%blk(b)%face(f)%eps * sigma * grid%blk(b)%face(f)%T**4
              if (Onwall > 2) then
                do j = 1, Ny; do k = 1, Nz
                  IOwall%block(bb)%vars(3:7, 1, j, k) = grid%blk(b)%face(f)%a(j, k, 0:Ngg)
                end do; end do
              end if

            case (3, 4)
              allocate(IOwall%block(bb)%mesh(1:3, 0:Nx, 0:0, 0:Nz))
              allocate(IOwall%block(bb)%vars(1:Onwall, 1:Nx, 1, 1:Nz))
              do k = 0, Nz; do i = 0, Nx
                IOwall%block(bb)%mesh(1, i, 0, k) = grid%blk(b)%x(i+1, merge(1,Ny+1,f==3), k+1)
                IOwall%block(bb)%mesh(2, i, 0, k) = grid%blk(b)%y(i+1, merge(1,Ny+1,f==3), k+1)
                IOwall%block(bb)%mesh(3, i, 0, k) = grid%blk(b)%z(i+1, merge(1,Ny+1,f==3), k+1)
              end do; end do
              IOwall%block(bb)%vars(1, 1:Nx, 1, 1:Nz) = &
                grid%blk(b)%face(f)%eps(1:Nx, 1:Nz)
              IOwall%block(bb)%vars(2, 1:Nx, 1, 1:Nz) = &
                grid%blk(b)%face(f)%Qtot * grid%blk(b)%face(f)%eps &
                - grid%blk(b)%face(f)%eps * sigma * grid%blk(b)%face(f)%T**4
              if (Onwall > 2) then
                do i = 1, Nx; do k = 1, Nz
                  IOwall%block(bb)%vars(3:7, i, 1, k) = grid%blk(b)%face(f)%a(i, k, 0:Ngg)
                end do; end do
              end if

            case (5, 6)
              allocate(IOwall%block(bb)%mesh(1:3, 0:Nx, 0:Ny, 0:0))
              allocate(IOwall%block(bb)%vars(1:Onwall, 1:Nx, 1:Ny, 1))
              do j = 0, Ny; do i = 0, Nx
                IOwall%block(bb)%mesh(1, i, j, 0) = grid%blk(b)%x(i+1, j+1, merge(1,Nz+1,f==5))
                IOwall%block(bb)%mesh(2, i, j, 0) = grid%blk(b)%y(i+1, j+1, merge(1,Nz+1,f==5))
                IOwall%block(bb)%mesh(3, i, j, 0) = grid%blk(b)%z(i+1, j+1, merge(1,Nz+1,f==5))
              end do; end do
              IOwall%block(bb)%vars(1, 1:Nx, 1:Ny, 1) = &
                grid%blk(b)%face(f)%eps(1:Nx, 1:Ny)
              IOwall%block(bb)%vars(2, 1:Nx, 1:Ny, 1) = &
                grid%blk(b)%face(f)%Qtot * grid%blk(b)%face(f)%eps &
                - grid%blk(b)%face(f)%eps * sigma * grid%blk(b)%face(f)%T**4
              if (Onwall > 2) then
                do i = 1, Nx; do j = 1, Ny
                  IOwall%block(bb)%vars(3:7, i, j, 1) = grid%blk(b)%face(f)%a(i, j, 0:Ngg)
                end do; end do
              end if
          end select
        end associate
      end do
    end do
  end subroutine Initialize_Wall_File


  subroutine Write_Wall_Solution(grid, file, fmt)
    use Lib_VTK
    use Lib_Tecplot
    implicit none
    type(GROOT_domain_type), intent(in)  :: grid
    character(len=*),        intent(in)  :: file
    character(len=*),        intent(in)  :: fmt
    character(llen) :: path, localpath_vtk
    integer         :: E_IO

    path = 'OUTPUT/'

    if (index(fmt, 'vtk') > 0) then
      if (index(fmt, 'ascii') > 0) then
        IOwall%vtk%format = 'ascii'
      else
        IOwall%vtk%format = 'raw'
      end if
      localpath_vtk = trim(path)//'vtk/'
      call execute_command_line('mkdir -p '//trim(localpath_vtk))
      E_IO = vtk_write_structured_multiblock(orion=IOwall, &
               vtspath=trim(localpath_vtk)//trim(file), &
               vtmpath=trim(path)//trim(file), &
               varnames=Owallnames, time=grid%time)
    else
      if (index(fmt, 'binary') > 0) then
        IOwall%tec%format = 'binary'
        IOwall%tec%extension = '.szplt'
      else
        IOwall%tec%format = 'ascii'
        IOwall%tec%extension = '.tec'
      end if
      E_IO = tec_write_structured_multiblock(orion=IOwall, varnames=Owallnames, &
               filename=trim(path)//trim(file)//'.tec')
    end if
  end subroutine Write_Wall_Solution

end module GROOT_IO_Wall
