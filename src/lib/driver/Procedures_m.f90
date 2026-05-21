module GROOT_Procedures_m
  use GROOT_Wrap_Setup
  use GROOT_Wrap_Solve
  use GROOT_Wrap_Postprocess

  implicit none

  type :: GROOT_type
  
  contains
    procedure, nopass :: setup       => GROOT_setup
    procedure, nopass :: solve       => GROOT_solve
    procedure, nopass :: postprocess => GROOT_postprocess
  end type GROOT_type

end module GROOT_Procedures_m
