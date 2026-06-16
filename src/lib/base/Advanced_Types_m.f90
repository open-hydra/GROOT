
module GROOT_Advanced_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use GROOT_Mod_Face,  only: SBface             ! structured-block face data
  use Lib_ORION_Data                            ! ORION_data

  implicit none

  type :: block_type
    integer                    :: dim(3)       !< Cell count per direction (ghosts excluded)
  end type block_type


  !! Gaseous phase object
  type :: obj_gas_phase
    real(R8), allocatable :: p(:,:,:)        !< Pressure
    real(R8), allocatable :: T(:,:,:)        !< Temperature
    real(R8), allocatable :: R(:,:,:)        !< Gas constant
    real(R8), allocatable :: rho(:,:,:,:)    !< Species partial densities
  end type obj_gas_phase

  !! ======================================================
  !! GROOT extensions: radiation-specific fields
  !! ======================================================

  !! Block extension – adds all GROOT-specific data on top of geometry
  type, extends(block_type) :: GROOT_block_type
    !!Carrier-phase primitive variables (kept as obj_gas_phase sub-type)
    type(obj_gas_phase)   :: gas_phase      !< p, T, R, rho(:) from CFD solver

    !!Face data (wall BC, heat flux, emissivity, WSGG weights)
    type(SBface) :: face(6)                 !< One SBface per block face

    !!Nodal coordinates, indexed 1:Nx+1 x 1:Ny+1 x 1:Nz+1
    !!      (node i,j,k in mesh notation = x(i+1,j+1,k+1))
    real(R8), allocatable :: x(:,:,:)       !< Node x-coordinates
    real(R8), allocatable :: y(:,:,:)       !< Node y-coordinates
    real(R8), allocatable :: z(:,:,:)       !< Node z-coordinates

    !! Radiative field
    real(R8), allocatable :: ka(:,:,:,:)    !< Absorption coefficient (0:Ngg)
    real(R8), allocatable :: a (:,:,:,:)    !< Spectral weights (0:Ngg)
    real(R8), allocatable :: kb(:,:,:,:)    !< Malkmus beta_eff (0:Ngg) — SNB full only
    real(R8), allocatable :: Ib(:,:,:)      !< Blackbody power
    real(R8), allocatable :: source(:,:,:)  !< Radiative source div(q)
    real(R8), allocatable :: G(:,:,:)       !< Irradiation
  end type GROOT_block_type

  !! ======================================================
  !! GROOT domain type 
  !! ======================================================
  type :: GROOT_domain_type
    real(R8)  :: time             !< Solution time
    integer   :: nb               !< Number of blocks
    integer   :: nbound           !< Number of boundary face-cells
    type(GROOT_block_type), dimension(:), allocatable :: blk  !< Block array
  end type GROOT_domain_type

  !! ======================================================
  !! Top-level simulation type (used by all driver modules)
  !!
  !! Mirrors GROOT_simulation_type pattern.
  !! ======================================================
  type :: GROOT_simulation_type
    type(GROOT_domain_type) :: domain  !< Computational domain (blocks + BC)
    type(ORION_data)        :: OCP     !< Carrier-phase I/O field
    type(ORION_data)        :: ORP     !< Radiation-phase I/O field
  end type GROOT_simulation_type

end module GROOT_Advanced_Types_m
