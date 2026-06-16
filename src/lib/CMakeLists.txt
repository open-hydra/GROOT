file(GLOB_RECURSE GROOT_SRC
    CONFIGURE_DEPENDS
    ${CMAKE_SOURCE_DIR}/src/lib/*.f90
)

# Bake the SNB data directory path into a Fortran module at configure time
configure_file(
    "${CMAKE_SOURCE_DIR}/src/lib/physics/snb_datadir.f90.in"
    "${CMAKE_CURRENT_BINARY_DIR}/snb_datadir.f90"
    @ONLY
)
list(APPEND GROOT_SRC "${CMAKE_CURRENT_BINARY_DIR}/snb_datadir.f90")

add_library(GROOTL ${GROOT_SRC})
add_library(GROOTL::GROOTL ALIAS GROOTL)
target_link_libraries(GROOTL PRIVATE FiNeR::FiNeR ORION)

# Pass the source-tree path as a compile-time macro (used by Mod_Parameters.f90)
target_compile_definitions(GROOTL PRIVATE
    "SRC_PATH=\"${CMAKE_SOURCE_DIR}/src/\"")

target_include_directories(GROOTL INTERFACE
    $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)

if(USE_OPENMP)
    set_target_properties(GROOTL PROPERTIES
        COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
        LINK_FLAGS    "${OpenMP_Fortran_FLAGS}")
endif()
