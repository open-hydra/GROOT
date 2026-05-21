function(add_fortran_executable target source)

  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE GROOTL FiNeR::FiNeR)

  target_include_directories(${target} INTERFACE
      $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

  if(USE_OPENMP)
    target_compile_options(${target} PRIVATE ${OpenMP_Fortran_FLAGS})
    target_link_options(${target}    PRIVATE ${OpenMP_Fortran_FLAGS})
  endif()

endfunction()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_fortran_executable(${solver} "main.f90")
add_fortran_executable(DocGen   "docgen.f90")
