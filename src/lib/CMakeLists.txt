file(GLOB_RECURSE GROOT_SRC
    CONFIGURE_DEPENDS
    ${CMAKE_SOURCE_DIR}/src/lib/*.f90
)

# Lib_slw.f90 (real) and Lib_slw_stub.f90 both define module GROOT_Lib_SLW.
# Remove both from the glob; exactly one is added back below depending on USE_RADLIB.
list(FILTER GROOT_SRC EXCLUDE REGEX "Lib_slw")

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

# SLW model: real wrapper + radlib C interface when enabled, otherwise the stub.
if(USE_RADLIB)
    target_sources(GROOTL PRIVATE
        "${CMAKE_SOURCE_DIR}/src/lib/physics/Lib_slw.f90"
        "${RADLIB_PATH}/src/fortran/c_interface.cc")
    # PUBLIC so the final Fortran executable inherits the C++ runtime deps.
    target_link_libraries(GROOTL PUBLIC radlib ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES})
    target_compile_definitions(GROOTL PRIVATE WITH_RADLIB)
else()
    target_sources(GROOTL PRIVATE
        "${CMAKE_SOURCE_DIR}/src/lib/physics/Lib_slw_stub.f90")
endif()

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
