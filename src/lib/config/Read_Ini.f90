module GROOT_Read_Ini
  use iso_fortran_env, only: I4 => int32, R8 => real64

  implicit none
  private
  public :: Read_Inifile

contains

  subroutine Read_Inifile(simulation)
    use GROOT_Advanced_Types_m,      only: GROOT_simulation_type
    use GROOT_Config_Types_m,        only: obj_dtm, obj_rad_model, obj_io
    use GROOT_Global_m
    use GROOT_Input_Registry
    use GROOT_Backend_INI,           only: Load_Ini
    use GROOT_Register_IO,           only: Register_IO
    use GROOT_Register_Discretization, only: Register_Discretization
    use GROOT_Register_Options,      only: Register_Options
    use GROOT_Register_Model,        only: Register_Model, Apply_Model_Species
    use GROOT_IO_Solution,           only: Input_Solution_Setup
    use GROOT_Input_Registry,        only: Validate_Registry
    use finer,                       only: file_ini
    implicit none
    type(GROOT_simulation_type), intent(inout) :: simulation

    type(file_ini)      :: fini
    character(len=64)   :: turb, soot
    integer             :: error
    character(len=1024) :: validation_out
    integer             :: i, i_CS
    !! ----------------------------------------------------------
    !! 1. Load input.ini
    call fini%load(filename='input.ini')

    !! ----------------------------------------------------------
    !! 2. Build registry and set defaults
    
    call Register_IO()
    call Register_Discretization()
    call Register_Options()
    call Register_Model()

    !! ----------------------------------------------------------
    !! 3. Populate registry-pointed variables from file
    call Load_Ini(fini)

    !! ----------------------------------------------------------
    !! 4. Validate registry (required params set, allowed values respected)
    validation_out = Validate_Registry()


    !! ----------------------------------------------------------
    !! 6. Setup I/O extension and background format
    call Input_Solution_Setup()

    !! ----------------------------------------------------------
    !! 7. Read molecular weights from phase.txt
    call Read_IdealGas_Properties(trim('INPUT/phase.txt'))

    !! ----------------------------------------------------------
    !! 8. Parse species list and validate
    call Apply_Model_Species()

    turb  = ''
    error = 0
    call fini%get(section_name='MOSE-Physics', option_name='turbulence', &
                  val=turb, error=error)
    nrans = 0
    if (index(trim(turb), 'SA')        > 0) nrans = 1
    if (index(trim(turb), 'Wilcox2006')> 0) nrans = 2
    if (index(trim(turb), 'SST')       > 0) nrans = 2
    if (index(trim(turb), 'SSGLRR')   > 0) nrans = 7

    call fini%get(section_name='MOSE-Physics', option_name='soot-generation', &
                  val=soot, error=error)
    
    nsoot = 0 
    if ((trim(soot)=='LL91').or.(trim(soot)=='LIN')) then 
      i_CS = -1 
      do i = 1, nsc
        if (cfd_species(i) == 'C(gr)') then
          i_CS = i
          exit
        end if
      end do
      if (i_CS == -1) then 
        nsoot = 2 
      else
        nsoot = 1
      end if
    endif

    !! ----------------------------------------------------------
    !! 10. Copy Nr and source flag into simulation config
    Nr     = obj_dtm%Nr
    source = obj_dtm%source

  end subroutine Read_Inifile


  !! ----------------------------------------------------------
  !! Read molecular weights and species names from phase.txt
  subroutine Read_IdealGas_Properties(wmfile)
    use GROOT_Global_m, only: nsc, wm_tab, cfd_species
    use strings,              only: parse
    implicit none
    character(len=*), intent(in) :: wmfile
    integer        :: ios, i, unitfile
    character(256) :: wholestring
    character(256) :: args(2)

    open(newunit=unitfile, file=trim(wmfile), status='old', iostat=ios)
    if (ios /= 0) return

    ios = 0;  nsc = -1
    read(unitfile, *)
    do while (ios == 0)
      read(unitfile, '(A)', iostat=ios) wholestring
      nsc = nsc + 1
    end do

    if (allocated(wm_tab))      deallocate(wm_tab)
    if (allocated(cfd_species)) deallocate(cfd_species)
    allocate(wm_tab(1:nsc))
    allocate(cfd_species(1:nsc))

    rewind(unitfile)
    read(unitfile, *)
    do i = 1, nsc
      read(unitfile, '(A)') wholestring
      call parse(wholestring, ' ', args)
      read(args(1), *) cfd_species(i)
      read(args(2), *) wm_tab(i)
    end do
    close(unitfile)

  end subroutine Read_IdealGas_Properties

end module GROOT_Read_Ini
