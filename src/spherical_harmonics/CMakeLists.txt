# Source for special functions
set(SRC_SPHHARM 
    ${CMAKE_CURRENT_SOURCE_DIR}/spherical_harmonics.f90)

set(SRC ${SRC} ${SRC_SPHHARM} PARENT_SCOPE)