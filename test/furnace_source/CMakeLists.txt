set(target furnace)

add_executable(${target} furnace.f90)
target_link_libraries(${target} PRIVATE GROOTL FiNeR::FiNeR)

target_include_directories(${target} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

if(USE_OPENMP)
  target_compile_options(${target} PRIVATE ${OpenMP_Fortran_FLAGS})
  target_link_options(${target}    PRIVATE ${OpenMP_Fortran_FLAGS})
endif()

set_target_properties(${target} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/test/furnace_source/bin)
