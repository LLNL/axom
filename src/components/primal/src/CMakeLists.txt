#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set( primal_headers

     ## geometry
     geometry/BoundingBox.hpp
     geometry/HyperSphere.hpp
     geometry/MortonIndex.hpp
     geometry/NumericArray.hpp
     geometry/Point.hpp
     geometry/Polygon.hpp
     geometry/Ray.hpp
     geometry/Segment.hpp
     geometry/Triangle.hpp
     geometry/Vector.hpp

     ## operators
     operators/clip.hpp
     operators/closest_point.hpp
     operators/intersect.hpp
     operators/orientation.hpp
     operators/squared_distance.hpp

     operators/detail/clip_impl.hpp
     operators/detail/intersect_impl.hpp

     ## spatial acceleration
     spatial_acceleration/BVHTree.hpp
     spatial_acceleration/UniformGrid.hpp
   )

## NOTE: dummy source since BLT does not support header only libraries
set( primal_sources
    primal.cpp
   )

blt_add_library(
    NAME                  primal
    HEADERS               ${primal_headers}
    SOURCES               ${primal_sources}
    HEADERS_OUTPUT_SUBDIR primal
    DEPENDS_ON            axom_utils slic
    )
