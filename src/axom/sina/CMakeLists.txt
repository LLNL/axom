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
# Note: Sina also optionally depends on Adiak when AXOM_USE_ADIAK=ON
#------------------------------------------------------------------------------
axom_component_requires(NAME       Sina
                        TPLS       Conduit )

#------------------------------------------------------------------------------
# Set sina version
#------------------------------------------------------------------------------
set(SINA_VERSION_MAJOR 1)
set(SINA_VERSION_MINOR 14)
set(SINA_VERSION_PATCH 0)
axom_configure_file ( config.hpp.in
                      ${PROJECT_BINARY_DIR}/include/axom/sina/config.hpp )

#------------------------------------------------------------------------------
# Specify the sina headers/sources
#------------------------------------------------------------------------------
set(sina_headers
    core/ConduitUtil.hpp
    core/Curve.hpp
    core/CurveSet.hpp
    core/DataHolder.hpp
    core/Datum.hpp
    core/Document.hpp
    core/File.hpp
    core/ID.hpp
    core/Record.hpp
    core/Relationship.hpp
    core/Run.hpp
    )

set(sina_sources
    core/ConduitUtil.cpp
    core/Curve.cpp
    core/CurveSet.cpp
    core/DataHolder.cpp
    core/Datum.cpp
    core/Document.cpp
    core/File.cpp
    core/ID.cpp
    core/Record.cpp
    core/Relationship.cpp
    core/Run.cpp
    )

# Add Adiak header and source
blt_list_append( TO sina_headers ELEMENTS core/AdiakWriter.hpp IF AXOM_USE_ADIAK )
blt_list_append( TO sina_sources ELEMENTS core/AdiakWriter.cpp IF AXOM_USE_ADIAK )

# Add fortran interface for Sina
if (ENABLE_FORTRAN)
  blt_list_append( TO sina_headers ELEMENTS interface/sina_fortran_interface.h)
  blt_list_append( TO sina_sources
                   ELEMENTS interface/sina_fortran_interface.cpp interface/sina_fortran_interface.f90)
endif()

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(sina_depends
    core
    conduit::conduit
    )

blt_list_append( TO sina_depends ELEMENTS adiak::adiak IF AXOM_USE_ADIAK )

axom_add_library(NAME       sina
                 SOURCES    ${sina_sources}
                 HEADERS    ${sina_headers}
                 DEPENDS_ON ${sina_depends}
                 FOLDER     axom/sina)

axom_write_unified_header(NAME    sina
                          HEADERS ${sina_headers})

axom_install_component(NAME    sina
                       HEADERS ${sina_headers})

#------------------------------------------------------------------------------
# Add tests and examples
#------------------------------------------------------------------------------
if(AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if(AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()
