# Dielectric averaging functions
set(SRC_DIELECTRIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/angular_average.f90)
set(SRC ${SRC} ${SRC_DIELECTRIC} PARENT_SCOPE)