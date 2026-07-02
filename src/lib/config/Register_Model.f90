module GROOT_Register_Model
  implicit none
  private
  public :: Register_Model, Apply_Model_Species

  !! Module-level string for the species list (e.g. "H2O,CO2").
  !! Parsed into rad_species(:) by Apply_Model_Species after Load_Ini.
  character(len=256), save, target, public :: model_species_str = 'none'

contains

  subroutine Register_Model()
    use GROOT_Input_Registry
    use GROOT_Config_Types_m, only: obj_rad_model
    use GROOT_Global_m
    implicit none
    character(len=:), allocatable :: sec
    sec = trim(codename)//'-Model'

    call reg%add(sec, 'model', obj_rad_model%model, 'const', &
        'Spectral model for gas-phase absorption coefficient', &
        'const , gray, wsgg , wsgg-H2O , wsgg-H2OCO2 , snbw , snb', .false.)

    call reg%add(sec, 'wall_emissivity', obj_rad_model%eps_wall, '1.0', &
        'Wall emissivity (uniform, applied to all opaque walls)', '>= 0', .true.)

    call reg%add(sec, 'k', obj_rad_model%k_user, '1.0', &
        'Absorption coefficient [m^-1] — used only when model = const', '> 0', .false.)

    call reg%add(sec, 'slw_ngray', obj_rad_model%slw_ngray, '4', &
        'Number of gray gases for the SLW model (model = slw)', '1,25', .false.)

    call reg%add(sec, 'p_ref', obj_rad_model%p_ref, '1.0e5', &
        'Reference pressure [Pa] for SLW wall emission', '> 0', .false.)

    call reg%add(sec, 'species', model_species_str, 'none', &
        'Comma-separated list of radiating species (e.g. "H2O,CO2")', '', .false.)

  end subroutine Register_Model


  !! Parse model_species_str → rad_species(:) and validate against
  !! the internal species database.  Must be called after Load_Ini.
  subroutine Apply_Model_Species()
    use GROOT_Config_Types_m, only: obj_rad_model
    use GROOT_Global_m, only: nscrad, nscdat, rad_species, ind_species
    implicit none
    integer :: is
    character(20) :: sc_database(4)
    integer :: js, check


    !! For 'const', 'snbw', 'snb', 'slw': species looked up internally, no INI species list needed
    if (trim(obj_rad_model%model) == 'const' .or. &
        trim(obj_rad_model%model) == 'snbw'  .or. &
        trim(obj_rad_model%model) == 'snb'   .or. &
        trim(obj_rad_model%model) == 'slw') then
      nscrad = 0
      if (allocated(rad_species)) deallocate(rad_species)
      allocate(rad_species(0))
      if (allocated(ind_species)) deallocate(ind_species)
      allocate(ind_species(0))
      return
    end if

    !! Count species from comma-separated string
    nscrad = 1
    do is = 1, len_trim(model_species_str)
      if (model_species_str(is:is) == ',') nscrad = nscrad + 1
    end do

    if (allocated(rad_species)) deallocate(rad_species)
    allocate(rad_species(1:nscrad))
    read(model_species_str, *) rad_species(1:nscrad)

    !! Species database
    sc_database(1) = 'H2O'
    sc_database(2) = 'CO2'
    sc_database(3) = 'CO'
    sc_database(4) = 'CH4'
    nscdat = size(sc_database)

    !! Validate each radiating species is in the database
    do is = 1, nscrad
      check = 0
      do js = 1, nscdat
        if (sc_database(js) == rad_species(is)) check = 1
      end do
      if (check == 0) then
        write(*, '(A,A,A)') ' ERROR: radiating species "', trim(rad_species(is)), &
            '" not found in database.'
        stop
      end if
    end do

    !! Build ind_species mapping
    if (allocated(ind_species)) deallocate(ind_species)
    allocate(ind_species(1:nscdat))
    ind_species(:) = 0
    do is = 1, nscdat
      check = 0
      do js = 1, nscrad
        if (sc_database(is) == rad_species(js)) check = 1
      end do
      if (check == 1) then
        !! Map to cfd_species index
        call find_cfd_index(sc_database(is), ind_species(is))
      end if
    end do

    !! Resolve WSGG alias: wsgg → wsgg-H2O or wsgg-H2OCO2
    if (trim(obj_rad_model%model) == 'wsgg') then
      obj_rad_model%model = 'wsgg-H2O'
      do is = 1, nscrad
        if (trim(rad_species(is)) == 'CO2') then
          obj_rad_model%model = 'wsgg-H2OCO2'
          exit
        end if
      end do
    end if

  end subroutine Apply_Model_Species


  subroutine find_cfd_index(name, idx)
    use GROOT_Global_m, only: nsc, cfd_species
    implicit none
    character(len=*), intent(in)  :: name
    integer,          intent(out) :: idx
    integer :: js
    idx = 0
    do js = 1, nsc
      if (trim(cfd_species(js)) == trim(name)) then
        idx = js
        return
      end if
    end do
  end subroutine find_cfd_index

end module GROOT_Register_Model
