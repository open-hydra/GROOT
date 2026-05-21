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

  # ORION
  if(NOT TARGET ORION)
    set(orionD "${ORION_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/ORION/")
    add_subdirectory(${orionD})
  endif()

  # FiNeR  (ini-file reader)
  if(NOT TARGET FiNeR)
    set(finerD "${FINER_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}/build/lib/third_party/FiNeR/")
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
endif()
