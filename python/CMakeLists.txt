# Find f2py
set(cmd1 "locate")
set(cmd1_arg "bin/f2py")
set(cmd2 "head")
set(cmd2_arg "-n 1")
execute_process(COMMAND ${cmd1} ${cmd1_arg} COMMAND ${cmd2} ${cmd2_arg} 
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  RESULT_VARIABLE f2py_result
  OUTPUT_VARIABLE f2py_exec OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "f2py found : ${f2py}")

string(REPLACE ";" " " objects_str "$<TARGET_OBJECTS:inverse_dielectric_objects>")
message(STATUS "OBJECT FILES: ${objects_str}")

if (GPU)
    set(LINKS ${spg_lib} inverse_dielectric ${MAGMA_DIR}/lib/libmagma.so)
else()
    set(LINKS ${spg_lib} inverse_dielectric lapack)
endif(GPU)

set(f2py_command "${f2py} -c ${objects_str} interface.f90 -m 
                  inverse_dielectric_wrapper -I ${CMAKE_Fortran_MODULE_DIRECTORY} --fcompiler=intelem")

message(STATUS "COMMAND F2PY: ${f2py_command}")

execute_process(
    COMMAND ${f2py_command}
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_STRIP_TRAILING_WHITESPACE
)