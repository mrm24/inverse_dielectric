# Library for linear algebra
if(GPU)
    set(SRC_LINALG ${CMAKE_CURRENT_SOURCE_DIR}/gpu_magma_t.f90)
else()
    set(SRC_LINALG ${CMAKE_CURRENT_SOURCE_DIR}/cpu_magma_t.f90)
endif(GPU)

set(SRC_LINALG ${SRC_LINALG}
    ${CMAKE_CURRENT_SOURCE_DIR}/linalg.f90)

set(SRC ${SRC} ${SRC_LINALG} PARENT_SCOPE)