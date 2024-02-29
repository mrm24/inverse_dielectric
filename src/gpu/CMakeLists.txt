# Library for device host register
if (GPU)
   set(SRC_GPU ${CMAKE_CURRENT_SOURCE_DIR}/gpu_register.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/gpu_register.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/gpu_unified_vendors_interface.f90
      ${CMAKE_CURRENT_SOURCE_DIR}/gpu_world_t_gpu.f90)
else()
   set(SRC_GPU ${CMAKE_CURRENT_SOURCE_DIR}/gpu_world_t_cpu.f90)
endif()
set(SRC ${SRC} ${SRC_GPU} PARENT_SCOPE)

