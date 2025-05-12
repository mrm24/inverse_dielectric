#Process the folders
# Constants
add_subdirectory(constants)

# Special functions
add_subdirectory(special_functions)

# Cell
add_subdirectory(crystal)

# Interpolation
add_subdirectory(interpolation)

# Spherical quadrature
add_subdirectory(quadrature)

# Linear algebra
add_subdirectory(linear_algebra)

# The IDieL object itself
add_subdirectory(idiel)

message(STATUS "FILES: ${SRC}")
set(MOD_DIR "${CMAKE_BINARY_DIR}/include")
file(MAKE_DIRECTORY ${MOD_DIR})
set(CMAKE_Fortran_MODULE_DIRECTORY ${MOD_DIR})
if(CMAKE_Fortran_COMPILER MATCHES "IntelLLVM")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -module ${MOD_DIR}")
endif()
add_library(IDieL_objects OBJECT ${SRC})
add_library(IDieL SHARED $<TARGET_OBJECTS:IDieL_objects>)
set_target_properties(IDieL PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)

# Add symbols to the static library so we only need to link against it and not its dependencies
if (NVIDIA OR AMD)
    target_link_libraries(IDieL PUBLIC ${MAGMA_DIR}/lib/libmagma.so)
endif()

set_target_properties(IDieL PROPERTIES
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

install(TARGETS IDieL DESTINATION ${INSTALL_DIR}/lib)
