!! ============================================================
!! GROOT DocGen — generates docs/user/input/input-parameters.md
!! from the live parameter registry.
!!
!! Usage:
!!   cd <repo-root>
!!   ./bin/DocGen
!!
!! Output: docs/user/input/input-parameters.md
!! ============================================================
program GROOT_docgen
  use GROOT_Input_Registry
  use GROOT_Register_IO,            only: Register_IO
  use GROOT_Register_Discretization,only: Register_Discretization
  use GROOT_Register_Options,       only: Register_Options
  use GROOT_Register_Model,         only: Register_Model
  use GROOT_Config_Types_m
  use GROOT_Global_m

  implicit none

  !! Minimal initialisation so registry add_* routines can write
  !! defaults into the config objects without accessing unallocated
  !! allocatable arrays.
  nsc  = 0
  nrans   = 0
  nsoot= 0
  if (.not. allocated(cfd_species)) allocate(cfd_species(0))
  if (.not. allocated(wm_tab))      allocate(wm_tab(0))
  if (.not. allocated(rad_species)) allocate(rad_species(0))
  if (.not. allocated(ind_species)) allocate(ind_species(0))

  !! Build the registry (sets defaults, registers all parameters)
  call Register_IO()
  call Register_Discretization()
  call Register_Options()
  call Register_Model()

  !! Write the Markdown table
  call reg%generate_markdown('docs/user/input/input-parameters.md')

  write(*,'(A)') ' DocGen: wrote docs/user/input/input-parameters.md'

end program GROOT_docgen
