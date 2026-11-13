module GROOT_IO_Solution
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Input_Solution_Setup
  public :: Output_Solution_Setup
  public :: read_wall_tec
  public :: vtk_varnames
  public :: IOtime

  integer          :: io_error
  real(R8)         :: IOtime = 0.0_R8   !! time of the gas field (the VTK files hold it; Tecplot ones do not)

  !! Procedure pointers (concrete target set by Input_Solution_Setup)
  procedure(r_solution_if), pointer, public :: read_bck      => null()
  procedure(w_solution_if), pointer, public :: write_solution => null()

  !! ------------------------------------------------------
  !! Abstract interfaces
  abstract interface
    subroutine r_solution_if(IOfield_gas)
      use Lib_ORION_Data
      implicit none
      type(orion_data), intent(inout)           :: IOfield_gas
    end subroutine r_solution_if

    subroutine w_solution_if(grid, IOfield, file)
      use Lib_ORION_Data
      use GROOT_Advanced_Types_m, only: GROOT_domain_type
      implicit none
      type(GROOT_domain_type), intent(in)    :: grid
      type(orion_data),        intent(inout) :: IOfield
      character(len=*),        intent(in)    :: file
    end subroutine w_solution_if
  end interface
  !! ------------------------------------------------------

contains

  !! @brief Set the input reader and resolve the MOSE output filename.
  subroutine Input_Solution_Setup()
    use GROOT_Config_Types_m, only: obj_io
    use GROOT_Parameters_m,   only: hlen
    use IR_Precision,         only: str
    implicit none
    logical         :: present
    integer         :: i
    character(len=hlen) :: try
    character(6)    :: extension

    if (index(obj_io%gas_format,'vtk')>0) then
      extension = '.vtm'
    else if (index(obj_io%gas_format,'ascii')>0) then
      extension = '.tec'
    else
      extension = '.szplt'
    end if
    read_bck => read_vtk_tec

    obj_io%nameinit = trim('OUTPUT/field')//trim(extension)
    inquire(file=obj_io%nameinit, exist=present)

    if (.not.present) then
      i = 0
      do
        i = i+1
        try = trim('OUTPUT/field')//trim(str(.true., i))//trim(extension)
        inquire(file=try,exist=present)
        if (present) then
          obj_io%nameinit = try
        else
          exit
        endif
      enddo
    endif

  end subroutine Input_Solution_Setup


  !! @brief Allocate IOfield vars and set procedure pointer for write_solution.
  subroutine Output_Solution_Setup(IOfield)
    use GROOT_Config_Types_m, only: obj_io
    use GROOT_Global_m,       only: model, Ngg
    use GROOT_Parameters_m,   only: hlen
    use IR_Precision,         only: str
    use Lib_ORION_Data,       only: orion_data
    implicit none
    type(orion_data), intent(inout) :: IOfield
    integer :: Onvar, b

    character(hlen) :: Varnames

    Varnames = ' '
    Onvar    = 0

    select case (trim(model))
      case ('gray', 'const', 'snbw', 'snb')
        Varnames = trim(Varnames)//'"Ib" "k" "divq"'
        Onvar    = 3
      case ('wsgg-H2O', 'wsgg-H2OCO2')
        Varnames = trim(Varnames)//'"Ib" "k(0)","k(1)","k(2)","k(3)","k(4)"'
        Varnames = trim(Varnames)//' "a(0)","a(1)","a(2)","a(3)","a(4)" "divq"'
        Onvar    = 12
    end select

    IOfield%vtk%node   = .false.
    IOfield%tec%node   = .false.


    !! VTK: 'ascii', or 'raw' (appended binary), as the wall file
    if (index(obj_io%sol_format,'ascii')>0) then
      IOfield%tec%extension = '.tec'
      IOfield%tec%format    = 'ascii'
      IOfield%vtk%format    = 'ascii'
    else
      IOfield%tec%extension = '.szplt'
      IOfield%tec%format    = 'binary'
      IOfield%vtk%format    = 'raw'
    endif

    IOfield%tec%bc     = .false.

    write_solution => write_vtk_tec

    if (Onvar /= size(IOfield%block(1)%vars, 1)) then
      do b = 1, size(IOfield%block)
        IOfield%block(b)%name = 'Block'//trim(str(.true., b))
        deallocate(IOfield%block(b)%vars)
        allocate(IOfield%block(b)%vars(1:Onvar, &
                                       1:IOfield%block(b)%Ni, &
                                       1:IOfield%block(b)%Nj, &
                                       1:IOfield%block(b)%Nk))
      end do
    end if
  end subroutine Output_Solution_Setup


  !! ----------------------------------------------------------
  subroutine write_vtk_tec(grid, IOfield, file)
    use GROOT_Config_Types_m,   only: obj_io
    use Lib_VTK
    use Lib_Tecplot
    use GROOT_Advanced_Types_m, only: GROOT_domain_type
    use GROOT_Global_m,         only: model, Ngg
    use GROOT_Parameters_m,     only: llen, hlen
    use IR_Precision,           only: str
    use Lib_ORION_Data,         only: orion_data
    implicit none
    type(GROOT_domain_type), intent(in)    :: grid
    type(orion_data),        intent(inout) :: IOfield
    character(len=*),        intent(in)    :: file
    ! Local
    character(llen) :: path, localpath_vtk
    character(hlen) :: Varnames
    integer :: b, i, E_IO

    Varnames = ' '
    path     = 'OUTPUT/'

    do b = 1, size(IOfield%block)
      IOfield%block(b)%name = 'Block'//trim(str(.true., b))

      select case (trim(model))
        case ('gray', 'const')
          associate(Nx => grid%blk(b)%dim(1), &
                    Ny => grid%blk(b)%dim(2), &
                    Nz => grid%blk(b)%dim(3))
            IOfield%block(b)%vars(1, 1:Nx, 1:Ny, 1:Nz) = grid%blk(b)%Ib   (1:Nx, 1:Ny, 1:Nz)
            IOfield%block(b)%vars(2, 1:Nx, 1:Ny, 1:Nz) = grid%blk(b)%ka   (1:Nx, 1:Ny, 1:Nz, 0)
            IOfield%block(b)%vars(3, 1:Nx, 1:Ny, 1:Nz) = grid%blk(b)%source(1:Nx, 1:Ny, 1:Nz)
          end associate
          Varnames = '"Ib" "k" "divq"'

        case ('snbw', 'snb')
          associate(Nx => grid%blk(b)%dim(1), &
                    Ny => grid%blk(b)%dim(2), &
                    Nz => grid%blk(b)%dim(3))
            IOfield%block(b)%vars(1, 1:Nx, 1:Ny, 1:Nz) = grid%blk(b)%Ib(1:Nx, 1:Ny, 1:Nz)
            IOfield%block(b)%vars(2, 1:Nx, 1:Ny, 1:Nz) = &
              sum(grid%blk(b)%ka(1:Nx, 1:Ny, 1:Nz, 1:Ngg) * &
                  grid%blk(b)%a (1:Nx, 1:Ny, 1:Nz, 1:Ngg), dim=4)
            IOfield%block(b)%vars(3, 1:Nx, 1:Ny, 1:Nz) = grid%blk(b)%source(1:Nx, 1:Ny, 1:Nz)
          end associate
          Varnames = '"Ib" "k" "divq"'

        case ('wsgg-H2O', 'wsgg-H2OCO2')
          associate(Nx => grid%blk(b)%dim(1), &
                    Ny => grid%blk(b)%dim(2), &
                    Nz => grid%blk(b)%dim(3))
            IOfield%block(b)%vars(1,  1:Nx, 1:Ny, 1:Nz) = grid%blk(b)%Ib(1:Nx, 1:Ny, 1:Nz)
            do i = 0, Ngg
              IOfield%block(b)%vars(2+i, 1:Nx, 1:Ny, 1:Nz) = grid%blk(b)%ka(1:Nx, 1:Ny, 1:Nz, i)
            end do
            do i = 0, Ngg
              IOfield%block(b)%vars(7+i, 1:Nx, 1:Ny, 1:Nz) = grid%blk(b)%a (1:Nx, 1:Ny, 1:Nz, i)
            end do
            IOfield%block(b)%vars(12, 1:Nx, 1:Ny, 1:Nz) = grid%blk(b)%source(1:Nx, 1:Ny, 1:Nz)
          end associate
          Varnames = '"Ib" "k(0)","k(1)","k(2)","k(3)","k(4)" "a(0)","a(1)","a(2)","a(3)","a(4)" "divq"'
      end select
    end do

    if (index(obj_io%sol_format,'vtk')>0) then
      localpath_vtk = trim(path)//'vtk/'
      call execute_command_line('mkdir -p '//trim(localpath_vtk))
      Varnames = vtk_varnames(Varnames)
      E_IO = vtk_write_structured_multiblock(orion=IOfield,vtspath=trim(localpath_vtk)//trim(file), &
                                                            vtmpath=trim(path)//trim(file),varnames=Varnames,time=grid%time)
    else
      E_IO = tec_write_structured_multiblock(orion=IOfield,varnames=Varnames,filename=trim(path)//trim(file)//trim(IOfield%tec%extension))
    end if
  end subroutine write_vtk_tec


  !! ----------------------------------------------------------
  subroutine read_vtk_tec(IOfield_gas)
    use GROOT_Config_Types_m, only: obj_io
    use GROOT_Parameters_m,   only: clen, hlen
    use Lib_ORION_Data,       only: orion_data
    use Lib_Tecplot
    use Lib_VTK
    use strings, only: parse
    implicit none
    type(orion_data), intent(inout)           :: IOfield_gas
    ! Local
    integer             :: error
    character(clen)     :: format(2)
    character(len=hlen) :: vtmpath, vtspath

    call parse(obj_io%gas_format,' ', format)

    select case(trim(format(1)))
    case('tecplot')
      IOfield_gas%tec%format = trim(format(2))
      error = tec_read_structured_multiblock(orion=IOfield_gas,filename=trim(obj_io%nameinit))
    case('vtk')
      IOfield_gas%vtk%format = trim(format(2))
      call resolve_vtk_input_paths(obj_io%nameinit, vtmpath, vtspath)
      error = vtk_read_structured_multiblock(orion=IOfield_gas, vtmpath=trim(vtmpath), vtspath=trim(vtspath), &
                                             time=obj_io%IOtime)
      IOtime = obj_io%IOtime
    end select

    if (error/=0) obj_io%error_message = "[ERROR] reading input file "//trim(obj_io%nameinit)

  end subroutine read_vtk_tec



  !! ----------------------------------------------------------
  !! The names of the variables as the VTK writer takes them: separated by one space, with no quotes. The Tecplot form,
  !! "Ib" "k(0)","k(1)", would put the quotes inside the name attribute of the file, which is not valid XML.
  function vtk_varnames(names) result(out)
    use GROOT_Parameters_m, only: hlen
    implicit none
    character(len=*), intent(in) :: names
    character(len=hlen)          :: out
    integer :: i, n

    out = ' '
    n   = 0
    do i = 1, len_trim(names)
      select case (names(i:i))
        case ('"')
          cycle
        case (',')
          n = n + 1
          out(n:n) = ' '
        case default
          n = n + 1
          out(n:n) = names(i:i)
      end select
    end do
    out = adjustl(out)
  end function vtk_varnames


  !! ----------------------------------------------------------
  !! Read wall.tec and locate the Tw variable index.
  !! For face zones (J=1 node) ORION sets block%Nj=0 and allocates vars
  !! with a zero-size j-dimension, so no var data is ever read.
  !! The second pass fixes this: it extracts zone names, computes nT,
  !! reallocates vars with the correct shape, and reads the Tw column.
  subroutine read_wall_tec(IOfield)
    use Lib_Tecplot
    use GROOT_Config_Types_m, only: obj_io
    use GROOT_Parameters_m,   only: hlen
    use Lib_ORION_Data,       only: orion_data
    use GROOT_Global_m,       only: nT
    implicit none
    type(orion_data), intent(inout) :: IOfield
    character(512)       :: line
    character(len=hlen)  :: wallfile
    integer        :: unitfree, i, j, ios, b
    integer        :: nvar_data       !! data vars = total vars - 3 coords
    integer        :: Ni_b, Nk_b     !! current block dimensions
    integer        :: idata           !! line counter within zone data section
    integer        :: skip_mesh       !! mesh lines before var data
    integer        :: skip_pre_Tw    !! var data lines before Tw column
    integer        :: Tw_lines        !! number of Tw values to read
    integer        :: k_idx, i_idx
    logical        :: in_data         !! .true. while reading a degenerate zone

    wallfile = 'OUTPUT/wall.tec'
    IOfield%tec%format = 'ascii'
    io_error = tec_read_structured_multiblock(orion=IOfield, filename=trim(wallfile))

    open(newunit=unitfree, file=trim(wallfile), status='old')
    !! nvar_data and nT from the variable names ORION has read, quoted or not (ORION 2 writes them bare)
    nT        = 0
    nvar_data = 0
    if (allocated(IOfield%varnames)) then
      nvar_data = size(IOfield%varnames) - 3
      nT        = findloc(IOfield%varnames, 'Tw', dim=1) - 3
    end if
    if (nT < 1) then
      write(*,'(A)') ' [ERROR] No variable Tw read from '//trim(wallfile)
      error stop
    end if
    b         = 0
    in_data   = .false.
    idata     = 0
    skip_mesh = 0;  skip_pre_Tw = 0;  Tw_lines = 0
    Ni_b = 0;  Nk_b = 0

    do
      read(unitfree, '(A)', iostat=ios) line
      if (ios /= 0) exit

      !! ZONE header — extract name; for Nj=0 blocks prepare manual Tw read
      if (index(line, 'ZONE') > 0 .and. index(line, 'ZONETYPE') == 0) then
        in_data = .false.
        b = b + 1
        if (b <= size(IOfield%block)) then
          i = index(line, 'T = ')
          if (i > 0) then
            i = i + 4
            j = index(line(i:), ',') - 1
            if (j <= 0) j = len_trim(line(i:))
            IOfield%block(b)%name = trim(adjustl(line(i:i+j-1)))
          end if
          !! ORION sets Nj=J-1=0 for face zones (J=1); vars is zero-size.
          !! Reallocate with the correct shape and read Tw manually.
          if (IOfield%block(b)%Nj == 0 .and. nvar_data > 0 .and. nT > 0) then
            Ni_b = IOfield%block(b)%Ni
            Nk_b = IOfield%block(b)%Nk           !! >= 1 (K=2 → Nk=1)
            if (allocated(IOfield%block(b)%vars)) deallocate(IOfield%block(b)%vars)
            allocate(IOfield%block(b)%vars(1:nvar_data, 1:Ni_b, 1:1, 1:Nk_b))
            IOfield%block(b)%vars = 0.0_R8
            IOfield%block(b)%Nj  = 1             !! j=1 access now valid
            !! ORION writes: ndir*(Ni+1)*1*(Nk+1) mesh lines then nvar_data*Ni*Nk var lines
            skip_mesh   = 3 * (Ni_b + 1) * (Nk_b + 1)
            skip_pre_Tw = (nT - 1) * Ni_b * Nk_b
            Tw_lines    = Ni_b * Nk_b
            idata       = 0
            in_data     = .true.
          end if
        end if
        cycle
      end if

      !! Manual data read for degenerate (Nj was 0) zones
      if (in_data) then
        idata = idata + 1
        if (idata > skip_mesh + skip_pre_Tw .and. &
            idata <= skip_mesh + skip_pre_Tw + Tw_lines) then
          j     = idata - skip_mesh - skip_pre_Tw
          k_idx = (j - 1) / Ni_b + 1
          i_idx = mod(j - 1, Ni_b) + 1
          read(line, *, iostat=ios) IOfield%block(b)%vars(nT, i_idx, 1, k_idx)
        end if
        if (idata >= skip_mesh + nvar_data * Ni_b * Nk_b) in_data = .false.
      end if
    end do
    close(unitfree)
  end subroutine read_wall_tec


  !! ----------------------------------------------------------
  !! The two paths of the VTK reader from the name of the .vtm file (OUTPUT/field.vtm): vtmpath is the file without its
  !! extension, vtspath the directory of the file (with its slash): the .vtm lists its block files, OUTPUT/vtk/fieldBlock1.vts,
  !! relative to its own directory
  subroutine resolve_vtk_input_paths(nameinit, vtmpath, vtspath)
    use GROOT_Parameters_m, only: hlen
    implicit none
    character(len=*),    intent(in)  :: nameinit
    character(len=hlen), intent(out) :: vtmpath, vtspath
    integer :: last_sep

    vtmpath  = nameinit(1:len_trim(nameinit)-4)
    last_sep = scan(trim(vtmpath), '/\', back=.true.)
    if (last_sep > 0) then
      vtspath = vtmpath(1:last_sep)
    else
      vtspath = ''
    end if
  end subroutine resolve_vtk_input_paths


end module GROOT_IO_Solution
