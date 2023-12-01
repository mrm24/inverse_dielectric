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
set_target_properties(inverse_dielectric_interfaces PROPERTIES
    INSTALL_RPATH "${MAGMA_DIR}/lib/libmagma.so"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

pybind11_add_module(InverseDielectric MODULE interface.cpp)
target_include_directories(InverseDielectric PRIVATE ${MOD_DIR} ${pybind11_INCLUDE_DIRS} ${Python3_INCLUDE_DIRS})
target_link_libraries(InverseDielectric PRIVATE inverse_dielectric_interfaces ${Python3_LIBRARIES})
set_target_properties(InverseDielectric PROPERTIES
    INSTALL_RPATH "${MAGMA_DIR}/lib/libmagma.so"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

# Install
install(TARGETS inverse_dielectric_interfaces DESTINATION ${INSTALL_DIR}/interfaces)
install(TARGETS InverseDielectric DESTINATION ${INSTALL_DIR}/interfaces)