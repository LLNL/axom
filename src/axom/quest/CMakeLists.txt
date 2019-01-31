#------------------------------------------------------------------------------
# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

################################
# Check necessary dependencies
################################
axom_component_requires(NAME       QUEST
                        COMPONENTS MINT PRIMAL SLAM SLIC)

################################
# Specify all headers/sources
################################
set( quest_headers

    SignedDistance.hpp

    ## All-nearest-neighbors query
    AllNearestNeighbors.hpp
    detail/AllNearestNeighbors_detail.hpp

    ## Geometry Headers
    geom/Brood.hpp
    geom/DenseOctreeLevel.hpp
    geom/ImplicitGrid.hpp
    geom/InOutOctree.hpp
    geom/OctreeBase.hpp
    geom/OctreeLevel.hpp
    geom/SparseOctreeLevel.hpp
    geom/SpatialOctree.hpp

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
     primal
     )

if (ENABLE_OPENMP)
    list(APPEND quest_depends_on openmp)
endif()

if (ENABLE_MPI)
    list(APPEND quest_headers stl/PSTLReader.hpp)
    list(APPEND quest_sources stl/PSTLReader.cpp)

    if(AXOM_ENABLE_LUMBERJACK)
        list(APPEND quest_depends_on lumberjack mpi)
    endif()
endif()

if (SHROUD_FOUND)
    add_subdirectory(interface)
endif()

# Add sparsehash when available
if(SPARSEHASH_FOUND)
    list(APPEND quest_depends_on sparsehash)
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

################################
# Make/Install the library
################################
blt_add_library(
    NAME        quest
    SOURCES     ${quest_sources}
    HEADERS     ${quest_headers}
    DEPENDS_ON  ${quest_depends_on}
    FOLDER      axom/quest
    )

axom_write_unified_header(NAME    quest
                          HEADERS ${quest_headers}
                          )

axom_install_component(NAME    quest
                       HEADERS ${quest_headers}
                       )

################################
# Add examples
################################
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# Add tests
################################
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# Add docs
################################
if (AXOM_ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Add code checks
################################
axom_add_code_checks(PREFIX quest)
