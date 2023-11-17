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
add_library(dielgamma SHARED ${SRC})


