#Process the folders

# Constants
add_subdirectory(constants)

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
add_library(inverse_dielectric SHARED $<TARGET_OBJECTS:inverse_dielectric_objects>)