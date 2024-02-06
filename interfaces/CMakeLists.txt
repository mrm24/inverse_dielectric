# This first creates a CXX library and then uses it to create a python version

find_package (Python3 COMPONENTS Interpreter Development)
include_directories(${Python3_INCLUDE_DIRS})

set(pybind11_DIR ${PYBIND_CMAKE_DIR})
find_package(pybind11 REQUIRED CONFIG)

set(SRC_INTERFACES interface.f90
                   interface.cpp)

set(MOD_DIR "${CMAKE_BINARY_DIR}/include")

add_library(IDieL_interfaces SHARED ${SRC_INTERFACES})
target_include_directories(IDieL_interfaces PUBLIC ${MOD_DIR} ${pybind11_INCLUDE_DIRS})
target_link_libraries(IDieL_interfaces IDieL ifcore)
set_target_properties(IDieL_interfaces PROPERTIES
    INSTALL_RPATH "${MAGMA_DIR}/lib/libmagma.so:${INSTALL_DIR}/lib/"
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

pybind11_add_module(IDieLPy MODULE interface.cpp)
target_include_directories(IDieLPy PRIVATE ${MOD_DIR} ${pybind11_INCLUDE_DIRS} ${Python3_INCLUDE_DIRS})
target_link_libraries(IDieLPy PRIVATE IDieL_interfaces ${Python3_LIBRARIES})
set_target_properties(IDieLPy PROPERTIES
    INSTALL_RPATH "${MAGMA_DIR}/lib/libmagma.so:${INSTALL_DIR}/lib:${INSTALL_DIR}/interfaces" 
    INSTALL_RPATH_USE_LINK_PATH TRUE
)

# Install
install(TARGETS IDieL_interfaces DESTINATION ${INSTALL_DIR}/interfaces)
install(TARGETS IDieLPy DESTINATION ${INSTALL_DIR}/interfaces)
