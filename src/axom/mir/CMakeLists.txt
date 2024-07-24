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
                        COMPONENTS SLIC SLAM PRIMAL
                        TPLS       Conduit)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(mir_headers
    MIRMesh.hpp
    MIRMeshTypes.hpp
    ZooClippingTables.hpp
    InterfaceReconstructor.hpp
    CellData.hpp
    MeshTester.hpp
    CellClipper.hpp
    MIRUtilities.hpp
    CellGenerator.hpp

    blueprint_utilities.hpp
    ClipField.hpp
    CoordsetBlender.hpp
    EquiZAlgorithm.hpp
    FieldBlender.hpp
    FieldSlicer.hpp
    MIRAlgorithm.hpp
    RecenterField.hpp
    utilities.hpp

    views/dispatch_coordset.hpp
    views/dispatch_material.hpp
    views/dispatch_rectilinear_topology.hpp
    views/dispatch_structured_topology.hpp
    views/dispatch_topology.hpp
    views/dispatch_uniform_topology.hpp
    views/dispatch_unstructured_topology.hpp
    views/dispatch_utilities.hpp
    views/ExplicitCoordsetView.hpp
    views/MaterialView.hpp
    views/NodeArrayView.hpp
    views/RectilinearCoordsetView.hpp
    views/Shapes.hpp
    views/StridedStructuredIndexing.hpp
    views/StructuredIndexing.hpp
    views/StructuredTopologyView.hpp
    views/UniformCoordsetView.hpp
    views/UnstructuredTopologyPolyhedralView.hpp
    views/UnstructuredTopologyMixedShapeView.hpp
    views/UnstructuredTopologySingleShapeView.hpp

    clipping/ClipCases.h
    clipping/ClipTableManager.hpp
    )

set(mir_sources
    MIRMesh.cpp
    InterfaceReconstructor.cpp
    ZooClippingTables.cpp
    CellData.cpp
    MeshTester.cpp
    CellClipper.cpp
    CellGenerator.cpp

    blueprint_utilities.cpp
    MIRAlgorithm.cpp
    EquiZAlgorithm.cpp
    clipping/ClipCasesHex.cpp
    clipping/ClipCasesPyr.cpp
    clipping/ClipCasesQua.cpp
    clipping/ClipCasesTet.cpp
    clipping/ClipCasesTri.cpp
    clipping/ClipCasesWdg.cpp
    views/MaterialView.cpp
    utilities.cpp
    )

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(mir_depends_on core slic slam primal conduit::conduit)

blt_add_library(
    NAME        mir
    SOURCES     ${mir_sources}
    HEADERS     ${mir_headers}
    DEPENDS_ON  ${mir_depends_on}
    FOLDER      axom/mir)

axom_write_unified_header(NAME    mir
                          HEADERS ${mir_headers} )

axom_install_component(NAME    mir
                       HEADERS ${mir_headers} )

#------------------------------------------------------------------------------
# Add tests, benchmarks and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()


#------------------------------------------------------------------------------
# Add code checks
#------------------------------------------------------------------------------
axom_add_code_checks(
        PREFIX mir )

