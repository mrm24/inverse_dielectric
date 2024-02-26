cmake_minimum_required(VERSION 3.20.0 FATAL_ERROR)

#Some options
option(GPU "GPU support" OFF)
option(BENCHMARK "Generate the benchmark executables" OFF)
option(INTERFACES "Python bindings" OFF)
option(SPGLIB "Compile with spglib" ON)
option(OMP5 "If OMP5 support is enabled it offload body corrections to the GPU : only gfortran >= 14.00 -- though it is slow" OFF)
option(DOUBLE "Use double precision in real and complex" ON)

# Set precision settings
if(NOT DOUBLE)
    message(STATUS "Compiling with single precision")
    add_compile_definitions(USE_SINGLE_PRECISION)
endif()

# Check GPU support
if(GPU)
    # This is only to automatically link GPU libraries
    project(inverse_dielectric_library Fortran CXX)
    add_compile_definitions(USE_GPU)
    if (OMP5)
      add_compile_definitions(HAVEOMP5)
    endif()
    # Setting MAGMA library
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
    set(MAGMA_F90_DIR "")
    set(SRC_MAGMA_F90 "")
endif()

# Check compilers for GPU offloading using OMP5
if(GPU AND OMP5)
 if (NOT CMAKE_Fortran_COMPILER MATCHES "gfortran" AND 
      CMAKE_Fortran_COMPILER_VERSION LESS 1400) 
   message(FATAL_ERROR "GPU offload of loops is only supported by gfortran version >= 14")
 endif()
endif()

# Select compilers
if(CMAKE_Fortran_COMPILER MATCHES "ifx")
    message(STATUS "Intel compilers")
    set(CMAKE_CXX_COMPILER "icpc")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -fPIC -qopenmp")
    set(CMAKE_Fortran_COMPILER "ifx")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -cpp -O3 -fPIC \
                            -no-wrap-margin -fp-model source -qopenmp \
                            ")
    # Set MKL as linear algebra library
    set(linalg_lib "mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m dl")
elseif(CMAKE_Fortran_COMPILER MATCHES "gfortran")
    message(STATUS "GNU compilers")
    set(CMAKE_CXX_COMPILER "g++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -fPIC -fopenmp")
    set(CMAKE_Fortran_COMPILER "gfortran")
    if(GPU AND OMP5)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -cpp -O3 -fPIC \
                                -ffree-line-length-512 -fopenmp \
                                -foffload-options=-lgfortran -foffload-options=-lm")
    else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -cpp -O3 -fPIC \
                                -ffree-line-length-512 -fopenmp \
                                ")
    endif()
    # Set OpenBLAS as linear algebra library
    set(linalg_lib "openblas")
elseif(CMAKE_Fortran_COMPILER MATCHES "nvfortran")
    message(STATUS "NVIDIA compilers")
    set(CMAKE_CXX_COMPILER "nvc++")
    set(CMAKE_Fortran_COMPILER "nvfortran")
    if (GPU)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -fPIC -mp=gpu")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -cpp -O3 -fPIC \
                              -mp=gpu \
                                ")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -fPIC -mp")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -cpp -O3 -fPIC \
                              -mp \
                                ")
    endif(GPU)
    set(linalg_lib "openblas")
else()
    message(FATAL_ERROR "Compiler is not recognized: only GNU, Intel and NVIDIA compilers are officially supported")
endif()

# Find SPGLIB
if (SPGLIB)
  add_compile_definitions(USE_SPGLIB)
  set(SPGLIBDIR "None" CACHE PATH "Path to spglib")
  set(CMAKE_PREFIX_PATH ${SPGLIBDIR})
  find_library(spg_lib NAMES spg spglib libsymspg symspg PATHS ${SPGLIBDIR} REQUIRED)
  message(STATUS "SPGLIB found : ${spg_lib}")
endif()

# Create INSTALL DIRECTORY
set(INSTALL_DIR "${PROJECT_SOURCE_DIR}/IDieL")
message(STATUS "INSTALL DIRECTORY : ${PROJECT_SOURCE_DIR}")
file(MAKE_DIRECTORY ${INSTALL_DIR})
file(MAKE_DIRECTORY ${INSTALL_DIR}/lib)

# Set the components source
set(SRC "${SRC_MAGMA_F90}")
set(MODDIR "")
add_subdirectory(src)

# Python bindings
if (INTERFACES)
  file(MAKE_DIRECTORY ${INSTALL_DIR}/interfaces)
  set(PYBIND_CMAKE_DIR "None" CACHE PATH "Path to CMake configuration of pybind")
  message(STATUS "Pybind11 CMake Configuration path : ${PYBIND_CMAKE_DIR}")
  add_subdirectory(interfaces)
endif(INTERFACES)

# Compile test
add_subdirectory(test)

# If selected compile benchmark programs
# it also downloads all needed data for the benchmark
if(BENCHMARK)
  add_subdirectory(benchmark)
endif(BENCHMARK)

# Documentation we generate the configuration with Doxygen if possible.
find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
  file(INSTALL ${CMAKE_CURRENT_SOURCE_DIR}/doc/logo.png DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
    COMMENT "Generating documentation with Doxygen" VERBATIM)
else()
  message(WARNING "Doxygen not found, documentation cannot be generated")
endif()
