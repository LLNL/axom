# Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sina -- API to write data to Sina's common file format
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify necessary dependencies
#
# Note: Sina also optionally depends on Adiak when AXOM_SINA_USE_ADIAK=ON
#------------------------------------------------------------------------------
axom_component_requires(NAME       Sina
                        TPLS       Conduit )

#------------------------------------------------------------------------------
# Add Sina options
#------------------------------------------------------------------------------
cmake_dependent_option( AXOM_SINA_USE_ADIAK "Enables Adiak support in Sina" ON
                       "AXOM_USE_ADIAK" OFF)

#------------------------------------------------------------------------------
# Generate Sina config
#------------------------------------------------------------------------------
axom_configure_file ( config.hpp.in ${PROJECT_BINARY_DIR}/include/axom/sina/config.hpp )

#------------------------------------------------------------------------------
# Specify the sina headers/sources
#------------------------------------------------------------------------------
set(sina_headers
    include/ConduitUtil.hpp
    include/CppBridge.hpp
    include/Curve.hpp
    include/CurveSet.hpp
    include/DataHolder.hpp
    include/Datum.hpp
    include/Document.hpp
    include/File.hpp
    include/ID.hpp
    include/Record.hpp
    include/Relationship.hpp
    include/Run.hpp
    )

set(sina_sources
    src/ConduitUtil.cpp
    src/Curve.cpp
    src/CurveSet.cpp
    src/DataHolder.cpp
    src/Datum.cpp
    src/Document.cpp
    src/File.cpp
    src/ID.cpp
    src/Record.cpp
    src/Relationship.cpp
    src/Run.cpp
    )

# Add Adiak header and source
blt_list_append( TO sina_headers ELEMENTS include/AdiakWriter.hpp IF AXOM_SINA_USE_ADIAK )
blt_list_append( TO sina_sources ELEMENTS src/AdiakWriter.cpp IF AXOM_SINA_USE_ADIAK )

# Add fortran interface for Sina
if (ENABLE_FORTRAN)
  blt_list_append( TO sina_headers ELEMENTS interface/include/sina_fortran_interface.h)
  blt_list_append( TO sina_sources
                   ELEMENTS interface/src/sina_fortran_interface.cpp interface/src/sina_fortran_interface.f90)
endif()

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(sina_depends
    core
    conduit::conduit
    )

blt_list_append( TO sina_depends ELEMENTS adiak::adiak IF AXOM_SINA_USE_ADIAK )

axom_add_library(NAME       sina
                 SOURCES    ${sina_sources}
                 HEADERS    ${sina_headers}
                 DEPENDS_ON ${sina_depends}
                 FOLDER     axom/sina)

axom_write_unified_header(NAME    sina
                          HEADERS ${sina_headers})

axom_install_component(NAME    sina
                       HEADERS ${sina_headers})

install(FILES           ${PROJECT_BINARY_DIR}/include/axom/sina/config.hpp
        DESTINATION     include/axom/sina
        )

#------------------------------------------------------------------------------
# Add tests and examples
#------------------------------------------------------------------------------
if(AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if(AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

#------------------------------------------------------------------------------
# Add code checks
#------------------------------------------------------------------------------
axom_add_code_checks(PREFIX sina)
