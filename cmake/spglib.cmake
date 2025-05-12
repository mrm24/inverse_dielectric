# This CMake file finds spglib, and in case it 
# is not found it decays to the prepacked version 

option(SPGLIB "Compile with spglib support" OFF)

if (SPGLIB)
  add_compile_definitions(USE_SPGLIB)
  find_library(libsymspg NAMES spg spglib libsymspg symspg PATHS ${SPGLIBDIR})
  if (NOT libsymspg)
    set(spglibInstallDir "${CMAKE_BINARY_DIR}/INTERNAL_spglib_install")
    file(MAKE_DIRECTORY ${spglibInstallDir})
    ExternalProject_Add(INTERNAL_SPGLIB
    SOURCE_DIR    "${CMAKE_SOURCE_DIR}/external/spglib-2.6.0/"
    BINARY_DIR    "${CMAKE_BINARY_DIR}/INTERNAL_spglib_build"
    CONFIGURE_COMMAND CC=${CMAKE_C_COMPILER}
                      CFLAGS=-fPIC FCCPP=cpp
                      ${CMAKE_COMMAND} -H<SOURCE_DIR> -B<BINARY_DIR>
                      -DSPGLIB_INSTALL=${spglibInstallDir} -DSPGLIB_WITH_TESTS=OFF
    BUILD_COMMAND   $(MAKE) clean && $(MAKE)
    INSTALL_COMMAND $(MAKE) install
    )

    link_libraries(libsymspg INTERFACE)

    add_library(libsymspg INTERFACE)
    target_link_libraries(libsymspg INTERFACE ${spglibInstallDir}/lib/libsymspg.a)
    add_dependencies(libsymspg INTERNAL_SPGLIB)
  endif()

  # Link the project against libsymspg
  link_libraries(libsymspg)

endif()
