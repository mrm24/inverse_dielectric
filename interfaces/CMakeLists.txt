# This first creates a CXX library and then uses it to create a python version

find_package (Python3 COMPONENTS Interpreter Development)
include_directories(${Python3_INCLUDE_DIRS})

set(pybind11_DIR ${PYBIND_CMAKE_DIR})
find_package(pybind11 REQUIRED CONFIG)

set(SRC_INTERFACES interface.f90
                   interface.cpp)

set(MOD_DIR "${CMAKE_BINARY_DIR}/include")

add_library(inverse_dielectric_interfaces SHARED ${SRC_INTERFACES})
target_include_directories(inverse_dielectric_interfaces PUBLIC ${MOD_DIR} ${pybind11_INCLUDE_DIRS})
target_link_libraries(inverse_dielectric_interfaces inverse_dielectric ifcore)

pybind11_add_module(InverseDielectric MODULE interface.cpp)
target_include_directories(InverseDielectric PRIVATE ${MOD_DIR} ${pybind11_INCLUDE_DIRS} ${Python3_INCLUDE_DIRS})
target_link_libraries(InverseDielectric PRIVATE inverse_dielectric_interfaces ${Python3_LIBRARIES})



