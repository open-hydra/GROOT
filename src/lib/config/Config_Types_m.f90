module GROOT_Config_Types_m
  use iso_fortran_env, only: I4 => int32, R8 => real64
  use GROOT_Parameters_m

  implicit none
  private

  !! ======================================================
  !! Simulation parameters
  !! ======================================================
  type :: simulation_parameters_t
    character(len=llen) :: warning_message 
    character(len=llen) :: error_message   
    character(len=llen) :: description    
    ! USER-DEFINED INPUTS
    integer  :: iter_threshold                    !< Max iterations
    ! Coupling flags
    logical  :: HYDRA_postprocess = .false.
    logical  :: wallcoupled       = .false.
    logical  :: fieldcoupled      = .false.
    ! Useful variables
    integer         :: nthreads         ! Number of threads for simulation
    real(R8)        :: cputime(2)       ! Simulation time duration
  end type simulation_parameters_t
  !! ------------------------------------------------------

  !! ======================================================
  !! DTM discretization options
  !! ======================================================
  type :: dtm_t
    character(len=llen) :: warning_message = ' '
    character(len=llen) :: error_message   = ' '
    character(len=llen) :: description     = ' '
    ! USER-DEFINED INPUTS
    integer  :: Nr                         !< Number of rays
    real(R8) :: tol                        !< Heat-flux convergence tolerance
    integer  :: itermax                    !< Max DTM iterations
    real(R8) :: wedge_angle_deg            !< Wedge angle (2-D axisym, deg)
    logical  :: source                     !< Compute radiative source term
    logical  :: twoDax                     !< 2-D axisymmetric flag
    logical  :: optically_thin             !< Optically thin approximation
  end type dtm_t
  !! ------------------------------------------------------

  !! ======================================================
  !! Radiation spectral model
  !! ======================================================
  type :: rad_model_t
    character(len=llen) :: warning_message = ' '
    character(len=llen) :: error_message   = ' '
    character(len=llen) :: description     = ' '
    ! USER-DEFINED INPUTS
    character(len=llen) :: model              !< Spectral model name
    real(R8)            :: eps_wall           !< Wall emissivity
    real(R8)            :: k_user             !< User-given abs. coeff (model='const')
  end type rad_model_t
  !! ------------------------------------------------------

  !! ======================================================
  !! Input / Output
  !! ======================================================
  type :: io_t
    character(len=llen) :: warning_message = ' '
    character(len=llen) :: error_message   = ' '
    character(len=llen) :: description     = ' '
    ! USER-DEFINED INPUTS
    character(len=clen) :: gas_format   !< Mesh format + encoding
    character(len=clen) :: sol_format   !< Output format + encoding
    ! Useful variables
    character(len=hlen) :: nameinit        !< Resolved gas-phase input filename
    character(4)        :: extension       !< File extension (.tec, .vtm, …)
    real(R8)            :: IOtime
  end type io_t
  !! ------------------------------------------------------

  public :: simulation_parameters_t, dtm_t, rad_model_t, io_t

  !! ======================================================
  !! Module-level singleton objects (accessed by name)
  !! ======================================================
  type(simulation_parameters_t),  public :: obj_sim_param
  type(dtm_t),        public :: obj_dtm
  type(rad_model_t),  public :: obj_rad_model
  type(io_t),         public :: obj_io
  !! ------------------------------------------------------

end module GROOT_Config_Types_m
