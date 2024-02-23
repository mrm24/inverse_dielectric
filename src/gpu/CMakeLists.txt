# Library for device host register
set(SRC_GPU ${CMAKE_CURRENT_SOURCE_DIR}/gpu_register.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/gpu_register.f90)
set(SRC ${SRC} ${SRC_GPU} PARENT_SCOPE)

