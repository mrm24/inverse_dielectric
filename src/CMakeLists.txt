#Process the folders

# Constants
add_subdirectory(constants)

# Special functions
add_subdirectory(spherical_harmonics)

# Cell
add_subdirectory(crystal)

# Spherical quadrature
add_subdirectory(quadrature)

# Linear algebra
add_subdirectory(linear_algebra)

# The dielectric average itself
add_subdirectory(dielectric)

message(STATUS "FILES: ${SRC}")
set(MOD_DIR "${CMAKE_BINARY_DIR}/include")
file(MAKE_DIRECTORY ${MOD_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${MOD_DIR})
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -module ${MOD_DIR}")
message(STATUS "MODULES dir : ${MOD_DIR}")
add_library(inverse_dielectric_objects OBJECT ${SRC})
add_library(inverse_dielectric STATIC $<TARGET_OBJECTS:inverse_dielectric_objects>)

# Add symbols to the static library so we only need to link against it and not its dependencies
if (GPU)
    set(LINKS ${spg_lib} mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl ${MAGMA_DIR}/lib/libmagma.so)
else()
    set(LINKS ${spg_lib} mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl)
endif(GPU)

target_link_libraries(inverse_dielectric ${LINKS})
set_target_properties(inverse_dielectric PROPERTIES
    INSTALL_RPATH "${MAGMA_DIR}/lib/libmagma.so"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

install(DIRECTORY ${MOD_DIR}
        DESTINATION ${INSTALL_DIR}
        DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                            GROUP_READ GROUP_EXECUTE
                            WORLD_READ WORLD_EXECUTE
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE
                         GROUP_READ
                         WORLD_READ
        )

install(TARGETS inverse_dielectric DESTINATION ${INSTALL_DIR}/lib)
