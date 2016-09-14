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
    geom/BVHTree.hpp
    geom/BoundingBox.hpp
    geom/Determinants.hpp
    geom/FullGridOctreeLevel.hpp
    geom/GridPointOctreeLevel.hpp
    geom/HyperSphere.hpp
    geom/InOutOctree.hpp
    geom/Intersection.hpp
    geom/MortonIndex.hpp
    geom/MortonOctreeLevel.hpp
    geom/NumericArray.hpp
    geom/OctreeBase.hpp
    geom/OctreeLevel.hpp
    geom/Orientation.hpp
    geom/Point.hpp
    geom/Ray.hpp
    geom/Segment.hpp
    geom/SpatialOctree.hpp
    geom/SquaredDistance.hpp
    geom/Triangle.hpp
    geom/Vector.hpp
    geom/closest_point.hpp
    
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
     common
     slic
     mint
     slam
     )

if (ENABLE_BOOST)
    list(APPEND quest_depends_on boost)
endif()

if (ENABLE_MPI)
  list(APPEND quest_headers stl/PSTLReader.hpp)  
  list(APPEND quest_sources stl/PSTLReader.cpp)
  list(APPEND quest_depends_on lumberjack)
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
    )
