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
    geom/MortonIndex.hpp
    geom/OctreeBase.hpp
    geom/OctreeLevel.hpp
    geom/SparseOctreeLevel.hpp
    geom/SpatialOctree.hpp

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
     )

if (BOOST_FOUND)
    list(APPEND quest_depends_on boost)
endif()

if (ENABLE_MPI)
  list(APPEND quest_headers stl/PSTLReader.hpp)
  list(APPEND quest_sources stl/PSTLReader.cpp)

  if(ENABLE_LUMBERJACK)
      list(APPEND quest_depends_on lumberjack)
  endif()
endif()

if(ENABLE_SHROUD)
    add_subdirectory(interface)
endif()

if (ENABLE_FORTRAN AND ENABLE_MPI)
  ## append shroud generated headers/sources to corresponding lists
  list(APPEND quest_headers interface/c_fortran/wrapQUEST.h)
  list(APPEND quest_sources interface/c_fortran/wrapQUEST.cpp)
  list(APPEND quest_sources interface/c_fortran/shroudrt.cpp)
  list(APPEND quest_sources interface/c_fortran/wrapfquest.f)
endif()

blt_add_library(
    NAME        quest
    SOURCES     ${quest_sources}
    HEADERS     ${quest_headers}
    HEADERS_OUTPUT_SUBDIR quest
    DEPENDS_ON  ${quest_depends_on}
    USE_OPENMP  ${ENABLE_OPENMP}
    )
