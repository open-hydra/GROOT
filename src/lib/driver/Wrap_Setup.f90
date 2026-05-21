module GROOT_Wrap_Setup
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: GROOT_setup

contains

  subroutine GROOT_setup(simulation)
    use GROOT_Advanced_Types_m,   only: GROOT_simulation_type
    use GROOT_Config_Types_m,     only: obj_sim_param, obj_dtm, obj_rad_model, obj_io
    use GROOT_Global_m
    use GROOT_Read_Ini,           only: Read_Inifile
    use GROOT_Assign_Setup,       only: Assign_Setup
    use GROOT_IO_Solution,        only: read_bck, read_wall_tec, Output_Solution_Setup
    use GROOT_Mod_Allocate_Data,  only: Setup_Data_Structure
    use GROOT_Mod_Metrics,        only: Setup_Metrics, Wedge, Import_Face_Coords
    use GROOT_Mod_Bound,          only: Setup_Bound
    use GROOT_Mod_Phase,          only: Setup_Gas, Setup_Rad, Setup_Wall
    use GROOT_IO_Wall,            only: Initialize_Wall_File
    use Lib_ORION_Data,           only: copyORION, orion_data
    implicit none
    type(GROOT_simulation_type), intent(inout) :: simulation
    type(orion_data) :: IOfield_wall
    integer :: ib

    !! ----------------------------------------------------------
    call Print_Header()

    !! Read input.ini → populate obj_* config objects
    call Read_Inifile(simulation)

    !! Sync config objects → GROOT_Global_m globals (for physics files)
    call Assign_Setup()

    !! Print simulation info
    call Print_Shell_Info()

    !! Load gas-phase field (mesh + primitive variables)
    call read_bck(simulation%OCP)

    !! Allocate block data structures (uses mesh info from OCP)
    call Setup_Data_Structure(simulation%domain, simulation%OCP)

    !! Compute mesh metrics (reads nodes from OCP)
    call Setup_Metrics(simulation%domain, simulation%OCP)

    !! Read wall.tec and set up wall BC data
    call read_wall_tec(IOfield_wall)
    call Setup_Bound(simulation%domain)

    !! Import gas phase and compute radiative properties
    call Setup_Gas(simulation%domain, simulation%OCP)

    call Setup_Rad(simulation%domain, simulation%ORP)
    call Setup_Wall(simulation%domain, IOfield_wall)

    !! Copy mesh from gas to rad IOfield
    do ib = 1, simulation%domain%Nb
      simulation%ORP%block(ib)%mesh = simulation%OCP%block(ib)%mesh
    end do

    !! Apply axisymmetric wedge correction, then re-sync face coordinates
    if (twoDax) then
      call Wedge(simulation%domain, wedge_angle_deg * pi / 180.0_R8)
      call Import_Face_Coords(simulation%domain)
    end if

    !! Prepare output IOfield
    call Output_Solution_Setup(simulation%ORP)

    !! ----------------------------------------------------------
    call Print_Loading_Summary()
    call Check_Input()
    call Stop_Simulation()

    call cpu_time(obj_sim_param%cputime(1))

  contains

    subroutine Print_Header()
      write(*,*)
      write(*,'(A)') ' _____  _____   ____   ____  _____'
      write(*,'(A)') '|  __ \|  __ \ / __ \ / __ \|_   _|'
      write(*,'(A)') '| |  \/ |__) | |  | | |  | | | |'
      write(*,'(A)') '| | __ |  _  /| |  | | |  | | | |'
      write(*,'(A)') '| |_\ \| | \ \| |__| | |__| |_| |_'
      write(*,'(A)') ' \____/|_|  \_\\____/ \____/|_____|'
      write(*,*)
      write(*,'(A)') '  Grid-based RadiatiOn fOrtran tOolbox'
      write(*,'(A)') '  DTM — Discrete Transfer Method'
      write(*,*)
    end subroutine Print_Header


    subroutine Print_Shell_Info()
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Setup'
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' DTM discretization'
      write(*,'(A,T35,I0)')  '   Number of rays',       obj_dtm%Nr
      write(*,'(A,T35,E10.3)') '   Flux tolerance',     obj_dtm%tol
      write(*,'(A,T35,I0)')  '   Max DTM iterations',   obj_dtm%itermax
      if (obj_dtm%twoDax) &
        write(*,'(A,T35,F8.3)') '   Wedge angle (deg)',  obj_dtm%wedge_angle_deg
      write(*,*)
      write(*,'(A)') ' Radiation model'
      write(*,'(A,T35,A)') '   Model',   trim(obj_rad_model%model)
      write(*,'(A,T35,F8.4)') '   Wall emissivity', obj_rad_model%eps_wall
      write(*,'(A)') ' ========================================================================================='
      write(*,*)
    end subroutine Print_Shell_Info


    subroutine Print_Loading_Summary()
      write(*,'(A)') ' ========================================================================================='
      write(*,'(A)') ' Loading'
      write(*,'(A)') ' ========================================================================================='
      if (len_trim(obj_rad_model%error_message) > 0) then
        write(*,'(A,T35,A)') '   Radiation model', 'FAIL'
        write(*,'(4X,A)') trim(obj_rad_model%error_message)
      else
        write(*,'(A,T35,A)') '   Radiation model', 'OK'
      end if
      if (len_trim(obj_io%error_message) > 0) then
        write(*,'(A,T35,A)') '   IO', 'FAIL'
        write(*,'(4X,A)') trim(obj_io%error_message)
      else
        write(*,'(A,T35,A)') '   IO', 'OK'
      end if
      write(*,'(A)') ' ========================================================================================='
    end subroutine Print_Loading_Summary


    subroutine Check_Input()
      use GROOT_Input_Registry, only: Validate_Registry
      implicit none
      character(len=1024) :: out
      out = Validate_Registry()
      if (index(out, 'ERROR') > 0) then
        write(*,'(A,T35,A)') '   Input file', 'FAIL'
        write(*,'(4X,A)') trim(out)
        stop
      else
        write(*,'(A,T35,A)') '   Input file', 'OK'
      end if
    end subroutine Check_Input


    subroutine Stop_Simulation()
      logical :: has_error
      has_error = .false.
      if (index(obj_rad_model%error_message, 'ERROR') > 0) has_error = .true.
      if (index(obj_io%error_message,        'ERROR') > 0) has_error = .true.
      if (has_error) stop
    end subroutine Stop_Simulation

  end subroutine GROOT_setup

end module GROOT_Wrap_Setup
