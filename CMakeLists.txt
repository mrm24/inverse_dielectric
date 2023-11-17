cmake_minimum_required(VERSION 3.6.0 FATAL_ERROR)
project(dielectric_average Fortran)

#Fortran flags
set(CMAKE_Fortran_COMPILER "ifort")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O3 -fPIC \
                         -no-wrap-margin -fp-model source -qopenmp \
                        ")
# set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O3 -fPIC -fbounds-check -Wunused \
#                         -Wconversion -Wunderflow -Wdo-subscript -fopenmp -ffree-form \
#                         -std=f2018")

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

#Compiles the components
set(SRC "")
add_subdirectory(src)

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