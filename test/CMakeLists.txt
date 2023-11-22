# Build all the necessary executables

set(linkdir ${CMAKE_Fortran_MODULE_DIRECTORY})

if (GPU)
    set(LINKS ${spg_lib} inverse_dielectric ${MAGMA_DIR}/lib/libmagma.so)
else()
    set(LINKS ${spg_lib} inverse_dielectric lapack)
endif(GPU)

add_executable(test_crystal test_crystal.f90)
target_include_directories(test_crystal PUBLIC ${linkdir})
target_link_libraries(test_crystal ${LINKS})


add_executable(test_quadrature_lebedev test_quadrature_lebedev.f90)
target_include_directories(test_quadrature_lebedev PUBLIC ${linkdir})
target_link_libraries(test_quadrature_lebedev ${LINKS})

add_executable(test_dielectric_average_iso test_dielectric_average_iso.f90)
execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/test_resources/dielectric/Binv.dat.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/test_resources/dielectric/head.dat.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/test_resources/dielectric/wingL.dat.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
execute_process( COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/test_resources/dielectric/wingU.dat.tar.gz WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(test_dielectric_average_iso PUBLIC ${linkdir})
target_link_libraries(test_dielectric_average_iso ${LINKS})

add_executable(test_inverse test_inverse.f90)
target_include_directories(test_inverse PUBLIC ${linkdir})
target_link_libraries(test_inverse ${LINKS})


# Enable testing
include(CTest)
# add_test(NAME test_crystal COMMAND test_crystal)
# add_test(NAME test_quadrature_lebedev COMMAND test_quadrature_lebedev)
# add_test(NAME test_dielectric_average_iso COMMAND test_dielectric_average_iso)
add_test(NAME test_inverse COMMAND test_inverse)
