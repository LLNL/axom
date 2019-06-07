# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# MIR -- Material Interface Reconstruction
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       MIR
                        COMPONENTS SLIC)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(mir_headers
    MIRMesh.hpp
    MIRMeshTypes.hpp
    ZooBitMaps.hpp
    InterfaceReconstructor.hpp
    CellData.hpp
    )

set(mir_sources
    ../Axom.cpp
    MIRMesh.cpp
    InterfaceReconstructor.cpp
    ZooBitMaps.cpp
    CellData.cpp
    )

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(mir_depends_on core slic)

blt_add_library(
    NAME        mir
    SOURCES     ${mir_sources}
    HEADERS     ${mir_headers}
    DEPENDS_ON  ${mir_depends_on}
    FOLDER      axom/mir
    OBJECT      TRUE )

axom_write_unified_header(NAME    mir
                          HEADERS ${mir_headers} )

axom_install_component(NAME    mir
                       HEADERS ${mir_headers} )

#------------------------------------------------------------------------------
# Add tests, benchmarks and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
  if (ENABLE_BENCHMARKS)
  #  add_subdirectory(benchmarks)
  endif()
endif()

if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()


#------------------------------------------------------------------------------
# Add code checks
#------------------------------------------------------------------------------
axom_add_code_checks(
        PREFIX mir )

