!! ============================================================
!! VTK input and output of GROOT.
!!
!! For each of the two VTK formats of the input file (sol-format: "vtk ascii", "vtk raw"):
!!
!!   1. Gas field in, as MOSE writes it: a two-block field is written with ORION (OUTPUT/field.vtm and its block
!!      files in OUTPUT/vtk/), read by GROOT (Input_Solution_Setup and read_bck) and compared with what was written:
!!      mesh, variables, number of blocks and the time of the field.
!!   2. Radiative solution out: a domain with known Ib, k and source is written by GROOT (write_solution) and read back
!!      with ORION: variables, mesh and time.
!!   3. Wall solution out: the six faces of the domain are written by GROOT (Write_Wall_Solution) and read back.
!!
!! The comparison with Tecplot (the gas field read from "tecplot ascii") shows that the two readers return the same
!! data. The program stops with an error status if any check fails. No external data is needed.
!! ============================================================
program vtk_io_test
  use iso_fortran_env,         only: R8 => real64
  use GROOT_Config_Types_m,    only: obj_rad_model, obj_io, obj_sim_param
  use GROOT_Assign_Setup,      only: Assign_Setup
  use GROOT_Advanced_Types_m,  only: GROOT_domain_type
  use GROOT_Global_m
  use GROOT_Mod_Allocate_Data, only: Setup_Data_Structure
  use GROOT_Mod_Metrics,       only: Setup_Metrics
  use GROOT_IO_Solution,       only: Input_Solution_Setup, Output_Solution_Setup, read_bck, write_solution, IOtime
  use GROOT_IO_Wall,           only: Initialize_Wall_File, Write_Wall_Solution, viscous_flag, IOwall
  use Lib_ORION_Data
  use Lib_VTK
  use Lib_Tecplot

  implicit none

  integer, parameter :: nvar_gas = 8
  real(R8), parameter :: time_gas = 12.5_R8, time_rad = 3.25_R8
  integer :: nfail = 0
  character(len=10), parameter :: formats(2) = [character(len=10) :: 'vtk ascii', 'vtk raw']
  integer :: n

  obj_sim_param%nthreads = 1
  obj_rad_model%model    = 'const'
  obj_rad_model%eps_wall = 1.0_R8
  obj_rad_model%k_user   = 1.0_R8
  nsc = 1; nsoot = 0; nrans = 0; nscrad = 0; nscdat = 0
  allocate(cfd_species(1:nsc), wm_tab(1:nsc), rad_species(0), ind_species(0))
  cfd_species(1) = 'N2'
  wm_tab(1)      = 28.0_R8
  call Assign_Setup()

  call execute_command_line('mkdir -p OUTPUT/vtk', wait=.true.)

  do n = 1, size(formats)
    write(*,'(/,A)') ' ==== '//trim(formats(n))//' ===='
    call test_gas_input(trim(formats(n)))
    call test_outputs(trim(formats(n)))
  end do
  write(*,'(/,A)') ' ==== same data from Tecplot ===='
  call test_tecplot_agrees()

  write(*,*)
  if (nfail > 0) then
    write(*,'(A,I0,A)') ' vtk_io: ', nfail, ' check(s) FAILED'
    error stop 1
  end if
  write(*,'(A)') ' vtk_io: all checks passed'

contains

  !! ----------------------------------------------------------
  subroutine check(ok, what)
    logical,          intent(in) :: ok
    character(len=*), intent(in) :: what
    if (ok) then
      write(*,'(A)') '   ok    '//what
    else
      write(*,'(A)') '   FAIL  '//what
      nfail = nfail + 1
    end if
  end subroutine check

  !! Largest difference of two arrays, relative to the largest value of the first
  function rel_diff(a, b) result(d)
    real(R8), intent(in) :: a(:), b(:)
    real(R8) :: d
    if (size(a) /= size(b)) then
      d = huge(1.0_R8)
    else if (size(a) == 0) then
      d = 0.0_R8
    else
      d = maxval(abs(a - b)) / max(maxval(abs(a)), tiny(1.0_R8))
    end if
  end function rel_diff

  !! The format accepts a relative difference: raw files hold the doubles themselves, the ascii ones their decimal form
  function tol_of(fmt) result(t)
    character(len=*), intent(in) :: fmt
    real(R8) :: t
    if (index(fmt, 'raw') > 0) then
      t = 0.0_R8
    else
      t = 1.0e-6_R8
    end if
  end function tol_of

  function vtk_format_of(fmt) result(f)
    character(len=*), intent(in) :: fmt
    character(len=5) :: f
    if (index(fmt, 'ascii') > 0) then
      f = 'ascii'
    else
      f = 'raw'
    end if
  end function vtk_format_of

  !! ----------------------------------------------------------
  !! A block of nx x ny x nz cells: a slightly deformed mesh and variables of any size, with the value of a variable
  !! and a cell given by a function that no layout error could reproduce
  subroutine make_block(field, b, name, nx, ny, nz, nvar, shift)
    type(orion_data), intent(inout) :: field
    integer,          intent(in)    :: b
    character(len=*), intent(in)    :: name
    integer,          intent(in)    :: nx, ny, nz, nvar
    real(R8),         intent(in)    :: shift
    integer :: i, j, k, s

    field%block(b)%name = name
    field%block(b)%Ni = nx; field%block(b)%Nj = ny; field%block(b)%Nk = nz
    allocate(field%block(b)%mesh(1:3, 0:nx, 0:ny, 0:nz), field%block(b)%vars(1:nvar, 1:nx, 1:ny, 1:nz))
    do k = 0, nz; do j = 0, ny; do i = 0, nx
      field%block(b)%mesh(1,i,j,k) = shift + 0.25_R8*i + 0.01_R8*j
      field%block(b)%mesh(2,i,j,k) = 0.30_R8*j + 0.02_R8*sin(0.5_R8*i)
      field%block(b)%mesh(3,i,j,k) = 0.20_R8*k + 0.01_R8*i
    end do; end do; end do
    do k = 1, nz; do j = 1, ny; do i = 1, nx; do s = 1, nvar
      field%block(b)%vars(s,i,j,k) = 100.0_R8*s + sin(0.37_R8*i + 0.91_R8*j + 1.70_R8*k + 0.13_R8*s)
    end do; end do; end do; end do
  end subroutine make_block

  subroutine make_gas_field(field)
    type(orion_data), intent(out) :: field
    allocate(field%block(2))
    call make_block(field, 1, 'Block1', 4, 3, 2, nvar_gas, 0.0_R8)
    call make_block(field, 2, 'Block2', 3, 2, 2, nvar_gas, 2.0_R8)
  end subroutine make_gas_field

  !! ----------------------------------------------------------
  !! 1. The field of the gas as MOSE writes it, read by GROOT
  subroutine test_gas_input(fmt)
    character(len=*), intent(in) :: fmt
    type(orion_data) :: wr, rd
    character(len=200) :: varnames
    integer :: b, e

    call make_gas_field(wr)
    wr%vtk%format = vtk_format_of(fmt)
    wr%vtk%node   = .false.
    varnames = 'rho u v w p T gamma R'
    e = vtk_write_structured_multiblock(vtspath='OUTPUT/vtk/field', vtmpath='OUTPUT/field', orion=wr, &
                                        varnames=varnames, time=time_gas)
    call check(e == 0, 'the field is written with ORION')

    obj_io%gas_format = fmt
    obj_io%error_message = ' '
    call Input_Solution_Setup()
    call check(trim(obj_io%nameinit) == 'OUTPUT/field.vtm', 'the input file is OUTPUT/field.vtm')
    call read_bck(rd)
    call check(len_trim(obj_io%error_message) == 0, 'the read reports no error')
    if (len_trim(obj_io%error_message) > 0) then
      write(*,'(A)') '         '//trim(obj_io%error_message)
      return
    end if

    call check(allocated(rd%block), 'the blocks are read')
    if (.not. allocated(rd%block)) return
    call check(size(rd%block) == 2, 'two blocks')
    if (size(rd%block) /= 2) return
    do b = 1, 2
      call check(rd%block(b)%Ni == wr%block(b)%Ni .and. rd%block(b)%Nj == wr%block(b)%Nj .and. &
                 rd%block(b)%Nk == wr%block(b)%Nk, 'block sizes')
      call check(rel_diff(pack(wr%block(b)%mesh, .true.), pack(rd%block(b)%mesh, .true.)) <= tol_of(fmt), &
                 'mesh of the block')
      call check(all(shape(wr%block(b)%vars) == shape(rd%block(b)%vars)), 'shape of the variables')
      call check(rel_diff(pack(wr%block(b)%vars, .true.), pack(rd%block(b)%vars, .true.)) <= tol_of(fmt), &
                 'variables of the block')
    end do
    call check(abs(IOtime - time_gas) <= 1.0e-12_R8, 'the time of the field is read')
  end subroutine test_gas_input

  !! ----------------------------------------------------------
  !! A domain of one block with known Ib, k, source, and its walls
  subroutine make_domain(grid, iogas)
    type(GROOT_domain_type), intent(out) :: grid
    type(orion_data),        intent(out) :: iogas
    integer :: i, j, k, ifa, nx, ny, nz

    allocate(iogas%block(1))
    call make_block(iogas, 1, 'Block1', 4, 3, 2, nvar_gas, 0.0_R8)
    call Setup_Data_Structure(grid, iogas)
    call Setup_Metrics(grid, iogas)

    nx = grid%blk(1)%dim(1); ny = grid%blk(1)%dim(2); nz = grid%blk(1)%dim(3)
    allocate(grid%blk(1)%ka    (1:nx, 1:ny, 1:nz, 0:Ngg))
    allocate(grid%blk(1)%a     (1:nx, 1:ny, 1:nz, 0:Ngg))
    allocate(grid%blk(1)%Ib    (1:nx, 1:ny, 1:nz))
    allocate(grid%blk(1)%source(1:nx, 1:ny, 1:nz))
    grid%blk(1)%a = 1.0_R8
    do k = 1, nz; do j = 1, ny; do i = 1, nx
      grid%blk(1)%Ib(i,j,k)      = 1.0e3_R8 + sin(0.37_R8*i + 0.91_R8*j + 1.7_R8*k)
      grid%blk(1)%ka(i,j,k,:)    = 0.5_R8 + 0.1_R8*cos(0.21_R8*i + 0.57_R8*j + 0.9_R8*k)
      grid%blk(1)%source(i,j,k)  = -2.0e4_R8 + 7.0_R8*sin(0.53_R8*i - 0.41_R8*j + 1.1_R8*k)
    end do; end do; end do
    grid%time = time_rad

    do ifa = 1, 6
      grid%blk(1)%face(ifa)%bc   = 301
      grid%blk(1)%face(ifa)%T    = 300.0_R8 + 10.0_R8*ifa
      grid%blk(1)%face(ifa)%eps  = 0.5_R8 + 0.05_R8*ifa
      grid%blk(1)%face(ifa)%Qtot = 1.0e4_R8 * ifa
      call set_face_pattern(grid, ifa)
    end do
    if (allocated(viscous_flag)) deallocate(viscous_flag)
    allocate(viscous_flag(1:1, 1:6))
    viscous_flag = .true.
  end subroutine make_domain

  !! The values of a face vary over its cells, so that a swapped index shows
  subroutine set_face_pattern(grid, ifa)
    type(GROOT_domain_type), intent(inout) :: grid
    integer,                 intent(in)    :: ifa
    integer :: p, q
    do q = 1, size(grid%blk(1)%face(ifa)%eps, 2)
      do p = 1, size(grid%blk(1)%face(ifa)%eps, 1)
        grid%blk(1)%face(ifa)%eps (p,q) = 0.5_R8 + 0.04_R8*sin(0.7_R8*p + 1.3_R8*q + ifa)
        grid%blk(1)%face(ifa)%Qtot(p,q) = 1.0e4_R8*ifa + 100.0_R8*sin(0.9_R8*p - 0.4_R8*q)
      end do
    end do
  end subroutine set_face_pattern

  !! ----------------------------------------------------------
  !! 2 and 3. The radiative solution and the walls, written by GROOT and read with ORION
  subroutine test_outputs(fmt)
    character(len=*), intent(in) :: fmt
    type(GROOT_domain_type) :: grid
    type(orion_data) :: iogas, iorad, rd
    real(R8) :: t
    integer  :: e, nx, ny, nz, bw, nface

    call make_domain(grid, iogas)
    nx = grid%blk(1)%dim(1); ny = grid%blk(1)%dim(2); nz = grid%blk(1)%dim(3)

    !! as Wrap_Setup prepares the radiative IOfield
    allocate(iorad%block(1))
    iorad%block(1)%Ni = nx; iorad%block(1)%Nj = ny; iorad%block(1)%Nk = nz
    allocate(iorad%block(1)%vars(1:4, 1:nx, 1:ny, 1:nz))
    iorad%block(1)%mesh = iogas%block(1)%mesh
    obj_io%sol_format = fmt
    call Output_Solution_Setup(iorad)
    call write_solution(grid, iorad, 'rad-field')

    rd%vtk%format = vtk_format_of(fmt)
    e = vtk_read_structured_multiblock(vtspath='OUTPUT/', vtmpath='OUTPUT/rad-field', orion=rd, time=t)
    call check(e == 0, 'the radiative solution is read back')
    if (e == 0) then
      call check(size(rd%block) == 1, 'one block')
      call check(abs(t - time_rad) <= 1.0e-12_R8, 'the time of the solution')
      call check(size(rd%block(1)%vars, 1) == 3, 'three variables: Ib, k, divq')
      call check(rel_diff(pack(grid%blk(1)%Ib, .true.), pack(rd%block(1)%vars(1,:,:,:), .true.)) <= tol_of(fmt), 'Ib')
      call check(rel_diff(pack(grid%blk(1)%ka(:,:,:,0), .true.), pack(rd%block(1)%vars(2,:,:,:), .true.)) &
                 <= tol_of(fmt), 'k')
      call check(rel_diff(pack(grid%blk(1)%source, .true.), pack(rd%block(1)%vars(3,:,:,:), .true.)) &
                 <= tol_of(fmt), 'divq')
      call check(rel_diff(pack(iogas%block(1)%mesh, .true.), pack(rd%block(1)%mesh, .true.)) <= tol_of(fmt), 'mesh')
    end if

    !! walls
    call Initialize_Wall_File(grid)
    call Write_Wall_Solution(grid, 'rad-wall', fmt)
    nface = size(IOwall%block)
    call check(nface == 6, 'six wall blocks')
    deallocate(rd%block)
    rd%vtk%format = vtk_format_of(fmt)
    e = vtk_read_structured_multiblock(vtspath='OUTPUT/', vtmpath='OUTPUT/rad-wall', orion=rd)
    call check(e == 0, 'the wall solution is read back')
    if (e == 0 .and. size(rd%block) == nface) then
      do bw = 1, nface
        !! ORION returns a face, a block with no cell in its normal direction, with its mesh and sizes but without its
        !! variables: they are read from the block file itself
        call check(rd%block(bw)%Ni == IOwall%block(bw)%Ni .and. rd%block(bw)%Nj == IOwall%block(bw)%Nj .and. &
                   rd%block(bw)%Nk == IOwall%block(bw)%Nk, 'wall block sizes of '//trim(IOwall%block(bw)%name))
        call check(rel_diff(pack(IOwall%block(bw)%mesh, .true.), pack(rd%block(bw)%mesh, .true.)) <= tol_of(fmt), &
                   'wall mesh of '//trim(IOwall%block(bw)%name))
        call check_wall_vars(fmt, bw)
      end do
    else if (e == 0) then
      call check(.false., 'as many wall blocks as written')
    end if
  end subroutine test_outputs

  !! ----------------------------------------------------------
  !! The values of the wall variables in the block file of a face, as written: the arrays after the Points, taken by their
  !! names (Owallnames of GROOT_IO_Wall: emissivity, Q), in ascii or in raw (appended) form
  subroutine check_wall_vars(fmt, bw)
    character(len=*), intent(in) :: fmt
    integer,          intent(in) :: bw
    character(len=*), parameter  :: names(2) = [character(len=10) :: 'emissivity', 'Q']
    character(len=:), allocatable :: buf, txt
    character(len=256) :: file
    real(R8), allocatable :: vals(:)
    integer :: u, fsize, s, n, i1, i2, p0, off, nbytes, ios, q

    file = 'OUTPUT/vtk/rad-wall'//trim(IOwall%block(bw)%name)//'.vts'
    open(newunit=u, file=trim(file), access='stream', form='unformatted', status='old', iostat=ios)
    call check(ios == 0, 'the block file '//trim(file)//' exists')
    if (ios /= 0) return
    inquire(unit=u, size=fsize)
    allocate(character(len=fsize) :: buf)
    read(u) buf
    close(u)
    !! txt: the same bytes with the line ends as blanks, for the text of the file (buf keeps the raw data as it is)
    txt = buf
    do i1 = 1, fsize
      if (txt(i1:i1) == achar(10) .or. txt(i1:i1) == achar(13)) txt(i1:i1) = ' '
    end do
    p0 = 0
    if (index(fmt, 'raw') > 0) then
      p0 = index(buf, 'encoding="raw">')
      p0 = p0 + index(buf(p0:), '_') - 1                  !! the byte before the first of the data
    end if

    do s = 1, 2
      n = size(IOwall%block(bw)%vars(s,:,:,:))
      allocate(vals(n))
      i1 = index(txt, 'Name="'//trim(names(s))//'"')
      call check(i1 > 0, 'array '//trim(names(s))//' of '//trim(IOwall%block(bw)%name))
      if (i1 > 0) then
        i2 = i1 + index(txt(i1:), '>') - 1                 !! end of its tag
        ios = 0
        if (index(fmt, 'raw') > 0) then
          off = i1 + index(txt(i1:i2), 'offset="') + 7      !! first digit of its offset
          q   = index(txt(off:i2), '"')
          read(txt(off:off+q-2), *, iostat=ios) off
          if (ios == 0) then
            nbytes = transfer(buf(p0+off+1:p0+off+4), nbytes)     !! the 4-byte size that precedes the data
            if (nbytes /= 8*n) ios = 1
          end if
          if (ios == 0) vals = transfer(buf(p0+off+5:p0+off+4+8*n), vals)
        else
          read(txt(i2+1:), *, iostat=ios) vals
        end if
        call check(ios == 0, 'values of '//trim(names(s))//' of '//trim(IOwall%block(bw)%name)//' are read')
        if (ios == 0) then
          call check(rel_diff(pack(IOwall%block(bw)%vars(s,:,:,:), .true.), vals) <= tol_of(fmt), &
                     'wall '//trim(names(s))//' of '//trim(IOwall%block(bw)%name))
        end if
      end if
      deallocate(vals)
    end do
  end subroutine check_wall_vars

  !! ----------------------------------------------------------
  !! The gas field, written in Tecplot ASCII and read by GROOT, is the one that the VTK file gives
  subroutine test_tecplot_agrees()
    type(orion_data) :: wr, rd
    character(len=200) :: varnames
    integer :: e, b

    call make_gas_field(wr)
    wr%tec%format = 'ascii'
    wr%tec%node   = .false.
    wr%tec%bc     = .false.
    varnames = '"rho" "u" "v" "w" "p" "T" "gamma" "R"'
    e = tec_write_structured_multiblock(orion=wr, varnames=varnames, filename='OUTPUT/field.tec')
    call check(e == 0, 'the field is written in Tecplot')
    obj_io%gas_format = 'tecplot ascii'
    obj_io%error_message = ' '
    call Input_Solution_Setup()
    call check(trim(obj_io%nameinit) == 'OUTPUT/field.tec', 'the input file is OUTPUT/field.tec')
    call read_bck(rd)
    call check(len_trim(obj_io%error_message) == 0 .and. allocated(rd%block), 'the Tecplot read reports no error')
    if (.not. allocated(rd%block)) return
    do b = 1, size(rd%block)
      call check(rel_diff(pack(wr%block(b)%vars, .true.), pack(rd%block(b)%vars, .true.)) <= 1.0e-6_R8, &
                 'Tecplot variables of the block')
    end do

    !! the wall file in Tecplot takes the extension of its format
    block
      type(GROOT_domain_type) :: grid
      type(orion_data)        :: iogas
      logical :: found
      call make_domain(grid, iogas)
      call Initialize_Wall_File(grid)
      call execute_command_line('rm -f OUTPUT/rad-wall.tec OUTPUT/rad-wall.szplt', wait=.true.)
      call Write_Wall_Solution(grid, 'rad-wall', 'tecplot ascii')
      inquire(file='OUTPUT/rad-wall.tec', exist=found)
      call check(found, 'the Tecplot ASCII wall file is OUTPUT/rad-wall.tec')
    end block
  end subroutine test_tecplot_agrees

end program vtk_io_test
