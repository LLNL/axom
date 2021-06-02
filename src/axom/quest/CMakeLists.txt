# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Quest component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       QUEST
                        COMPONENTS MINT PRIMAL SLAM SLIC)

#------------------------------------------------------------------------------
# Specify headers, sources and dependencies
#------------------------------------------------------------------------------
set( quest_headers

    SignedDistance.hpp

    ## All-nearest-neighbors query
    AllNearestNeighbors.hpp
    detail/AllNearestNeighbors_detail.hpp

    ## In/out query
    InOutOctree.hpp
    detail/inout/BlockData.hpp
    detail/inout/MeshWrapper.hpp
    detail/inout/InOutOctreeMeshDumper.hpp
    detail/inout/InOutOctreeStats.hpp
    detail/inout/InOutOctreeValidator.hpp

    # Mesh tester headers
    MeshTester.hpp

    # PointInCell headers
    PointInCell.hpp
    detail/PointFinder.hpp
    detail/PointInCellMeshWrapper_mfem.hpp

    ## STL library headers
    stl/STLReader.hpp

    ## quest interface headers
    interface/internal/mpicomm_wrapper.hpp
    interface/internal/QuestHelpers.hpp

    interface/inout.hpp
    interface/signed_distance.hpp

    )

set( quest_sources

    ## All-nearest-neighbors query
    AllNearestNeighbors.cpp

    ## Geometry Sources
    # <nothing yet>

    ## Mesh tester sources
    MeshTester.cpp

    ## STL library sources
    stl/STLReader.cpp

    ## quest interface sources
    interface/internal/QuestHelpers.cpp

    interface/inout.cpp
    interface/signed_distance.cpp

   )

set( quest_depends_on
     core
     slic
     mint
     slam
     spin
     primal
     fmt
     )

blt_list_append(TO quest_depends_on IF SPARSEHASH_FOUND ELEMENTS sparsehash)
blt_list_append(TO quest_depends_on IF ENABLE_OPENMP ELEMENTS openmp)
blt_list_append(TO quest_depends_on IF ENABLE_CUDA ELEMENTS cuda)

if (ENABLE_MPI)
    list(APPEND quest_headers stl/PSTLReader.hpp)
    list(APPEND quest_sources stl/PSTLReader.cpp)

    blt_list_append(TO quest_depends_on
                    IF AXOM_ENABLE_LUMBERJACK
                    ELEMENTS lumberjack )

    list(APPEND quest_depends_on mpi)
endif()

if (SHROUD_FOUND)
    add_subdirectory(interface)
endif()

if (ENABLE_FORTRAN)
    ## append shroud generated headers/sources to corresponding lists
    list(APPEND quest_headers interface/c_fortran/wrapQUEST.h)
    list(APPEND quest_headers interface/c_fortran/typesQUEST.h)
    list(APPEND quest_sources interface/c_fortran/wrapQUEST.cpp)
    list(APPEND quest_sources interface/c_fortran/wrapfquest.F)
endif()

if (BUILD_SHARED_LIBS)
    # Install setup.py file used to build quest Python module.
    configure_file(interface/python/setup.py.in
                   ${CMAKE_CURRENT_BINARY_DIR}/interface/python/setup.py)
    configure_file(interface/python/quest_test.py.in
                   ${CMAKE_CURRENT_BINARY_DIR}/interface/python/quest_test.py)
endif()

#------------------------------------------------------------------------------
# Make and install the library
#------------------------------------------------------------------------------
blt_add_library(
    NAME        quest
    SOURCES     ${quest_sources}
    HEADERS     ${quest_headers}
    DEPENDS_ON  ${quest_depends_on}
    FOLDER      axom/quest
    OBJECT      TRUE
    )

axom_write_unified_header(NAME    quest
                          HEADERS ${quest_headers}
                          )

axom_install_component(NAME    quest
                       HEADERS ${quest_headers}
                       )

#------------------------------------------------------------------------------
# Add tests and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

axom_add_code_checks(PREFIX quest)
