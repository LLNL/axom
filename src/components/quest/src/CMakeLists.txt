#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set( quest_headers

    SignedDistance.hpp

    ## Geometry Headers
    geom/Brood.hpp
    geom/DenseOctreeLevel.hpp
    geom/InOutOctree.hpp
    geom/OctreeBase.hpp
    geom/OctreeLevel.hpp
    geom/SparseOctreeLevel.hpp
    geom/SpatialOctree.hpp
    geom/MeshTester.hpp

    geom/detail/MeshTester_impl.hpp

    ## STL library headers
    stl/STLReader.hpp

    ## quest interface headers
    interface/quest.hpp

   )

set( quest_sources

    ## Geometry Sources
    # <nothing yet>

    ## STL library sources
    stl/STLReader.cpp

    ## quest interface sources
    interface/quest.cpp

   )

set( quest_depends_on
     axom_utils
     slic
     mint
     slam
     primal
     )

if (BOOST_FOUND)
    list(APPEND quest_depends_on boost)
endif()

if (ENABLE_OPENMP)
    list(APPEND quest_depends_on openmp)
endif()

if (ENABLE_MPI)
    list(APPEND quest_headers stl/PSTLReader.hpp)
    list(APPEND quest_sources stl/PSTLReader.cpp)

    if(ENABLE_LUMBERJACK)
        list(APPEND quest_depends_on lumberjack mpi)
    endif()
endif()

if (SHROUD_FOUND)
    add_subdirectory(interface)
endif()

if (ENABLE_FORTRAN)
    if (ENABLE_MPI)
        ## append shroud generated headers/sources to corresponding lists
        list(APPEND quest_headers interface/c_fortran/wrapQUEST.h)
        list(APPEND quest_sources interface/c_fortran/wrapQUEST.cpp)
        list(APPEND quest_sources interface/c_fortran/shroudrt.cpp)
        list(APPEND quest_sources interface/c_fortran/wrapfquest.f)
    endif()
endif()

blt_add_library(
    NAME        quest
    SOURCES     ${quest_sources}
    HEADERS     ${quest_headers}
    HEADERS_OUTPUT_SUBDIR quest
    DEPENDS_ON  ${quest_depends_on}
    )
install(TARGETS              quest
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include/quest
        )
install(FILES       ${quest_headers}
        DESTINATION include/quest
        )
if (ENABLE_FORTRAN)
    install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/axom_quest.mod DESTINATION lib/fortran)
endif()
