# CMake project file for 'GROOT'

cmake_minimum_required(VERSION 3.22)

# Name the code
set(solver GROOT)

project(${solver} LANGUAGES Fortran CXX C)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

  set(CMAKE_VERBOSE_MAKEFILE FALSE)
  # Search path for cmake modules
  set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

  # Parallelization options
  OPTION(USE_MPI    "Use the MPI library for parallelization" OFF)
  OPTION(USE_OPENMP "Use OpenMP for parallelization"          OFF)

  INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake)
  INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)

  # Set directories
  include(GNUInstallDirs)
  set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

  # Helper: initialize a git submodule from the correct git root
  macro(ensure_submodule_initialized SUBMOD_PATH)
    if (NOT EXISTS "${SUBMOD_PATH}/CMakeLists.txt")
      message(STATUS "Submodule at ${SUBMOD_PATH} not initialized. Attempting git submodule update --init...")
      get_filename_component(_sm_parent "${SUBMOD_PATH}" DIRECTORY)
      execute_process(
        COMMAND git rev-parse --show-toplevel
        WORKING_DIRECTORY "${_sm_parent}"
        OUTPUT_VARIABLE _sm_git_root
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE _sm_git_result)
      if(_sm_git_result EQUAL 0)
        file(RELATIVE_PATH _sm_rel "${_sm_git_root}" "${SUBMOD_PATH}")
        execute_process(
          COMMAND git submodule update --init "${_sm_rel}"
          WORKING_DIRECTORY "${_sm_git_root}"
          RESULT_VARIABLE _sm_init_result)
        if(NOT _sm_init_result EQUAL 0)
          message(FATAL_ERROR "Failed to initialize submodule: ${SUBMOD_PATH}")
        endif()
      else()
        message(FATAL_ERROR "Could not find git root for submodule: ${SUBMOD_PATH}")
      endif()
    endif()
  endmacro()

  # Path fallbacks (standalone builds)
  if (NOT DEFINED ORION_PATH)
    set(ORION_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/ORION/")
  endif()
  if (NOT DEFINED FINER_PATH)
    set(FINER_PATH "${CMAKE_CURRENT_SOURCE_DIR}/lib/third_party/FiNeR/")
  endif()

  # ORION
  if(NOT TARGET ORION)
    ensure_submodule_initialized("${ORION_PATH}")
    set(orionD "${ORION_PATH}" "${CMAKE_BINARY_DIR}/lib/ORION")
    add_subdirectory(${orionD})
  endif()

  # FiNeR  (ini-file reader)
  if(NOT TARGET FiNeR)
    ensure_submodule_initialized("${FINER_PATH}")
    set(finerD "${FINER_PATH}" "${CMAKE_BINARY_DIR}/lib/third_party/FiNeR")
    add_subdirectory(${finerD})
  endif()

endif()

# Build library
add_subdirectory(src/lib)

# Build application (only if this is the top-level project)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  add_subdirectory(src/app)
  add_subdirectory(test/homcube)
  add_subdirectory(test/cube_periodic)
  add_subdirectory(test/finite_cylinder)
  add_subdirectory(test/furnace_source)
  add_subdirectory(test/quadrilateral)
  add_subdirectory(test/emissivity)
  add_subdirectory(test/nozzle)
endif()
