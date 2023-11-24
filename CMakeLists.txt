cmake_minimum_required(VERSION 3.20.0 FATAL_ERROR)

#Some options
option(GPU "GPU support" OFF)
option(INTERFACES "Python bindings" OFF)

# Check GPU support
if(GPU)
    project(inverse_dielectric_library Fortran CUDA CXX)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED on)
    set(CMAKE_CXX_COMPILER "icpc")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC -qopenmp")
    set(CMAKE_Fortran_COMPILER "ifort")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -cpp -O3 -fPIC \
                            -no-wrap-margin -fp-model source -qopenmp \
                            ")
    add_compile_definitions(USE_GPU) 
    set(MAGMA_F90_DIR ${PROJECT_SOURCE_DIR}/external/magma)
    set(SRC_MAGMA_F90 
           ${MAGMA_F90_DIR}/magma2.F90
           ${MAGMA_F90_DIR}/magma2_common.F90
           ${MAGMA_F90_DIR}/magma2_sfortran.F90
           ${MAGMA_F90_DIR}/magma2_cfortran.F90
           ${MAGMA_F90_DIR}/magma2_dfortran.F90
           ${MAGMA_F90_DIR}/magma2_zfortran.F90)
    set(MAGMA_DIR "None" CACHE PATH "Path to magma library")
    message(STATUS "magma library path : ${MAGMA_DIR}")
else()
    project(inverse_dielectric_library Fortran CXX)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED on)
    set(CMAKE_CXX_COMPILER "icpc")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC -qopenmp")
    set(CMAKE_Fortran_COMPILER "ifort")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -cpp -O3 -fPIC \
                            -no-wrap-margin -fp-model source -qopenmp \
                            ")
    set(MAGMA_F90_DIR "")
    set(SRC_MAGMA_F90 "")
endif(GPU)

# Find SPGLIB
set(cmd1 "locate")
set(cmd1_arg "libsymspg.so")
set(cmd2 "head")
set(cmd2_arg "-n 1")
execute_process(COMMAND ${cmd1} ${cmd1_arg} COMMAND ${cmd2} ${cmd2_arg} 
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  RESULT_VARIABLE spg_result
  OUTPUT_VARIABLE spg_lib OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "SPGLIB found : ${spg_lib}")

#Set the components source
set(SRC "${SRC_MAGMA_F90}")
set(MODDIR "")
add_subdirectory(src)

# Python bindings
if (INTERFACES)
  set(PYBIND_CMAKE_DIR "None" CACHE PATH "Path to CMake configuration of pybind")
  message(STATUS "Pybind11 CMake Configuration path : ${MAGMA_DIR}")
  add_subdirectory(interfaces)
endif(INTERFACES)

# Compile test
add_subdirectory(test)

# Documentation we generate the configuration with Doxygen if possible.
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
  file(INSTALL ${CMAKE_CURRENT_SOURCE_DIR}/doc/greenx_logo.png DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
    COMMENT "Generating documentation with Doxygen" VERBATIM)
else()
  message(WARNING "Doxygen not found, documentation cannot be generated")
endif()
