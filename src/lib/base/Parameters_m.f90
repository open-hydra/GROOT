module GROOT_Parameters_m

  implicit none

  integer, parameter :: hlen = 1024
  integer, parameter :: llen = 256
  integer, parameter :: clen = 16

  character(len=clen) :: codename = 'GROOT'

  character(len=clen) :: GROOT_phase_prefix = 'rad-'
  character(len=clen) :: MOSE_phase_prefix  = ''

end module GROOT_Parameters_m
