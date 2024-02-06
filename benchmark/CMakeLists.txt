# Build all the necessary executables
set(MOD_DIR "${CMAKE_BINARY_DIR}/include")
set(linkdir ${MOD_DIR})

set(LINKS IDieL)

add_executable(procedures_3d procedures_3d.f90)
target_include_directories(procedures_3d PUBLIC ${linkdir})
target_link_libraries(procedures_3d ${LINKS})

add_executable(procedures_2d procedures_2d.f90)
target_include_directories(procedures_2d PUBLIC ${linkdir})
target_link_libraries(procedures_2d ${LINKS})

