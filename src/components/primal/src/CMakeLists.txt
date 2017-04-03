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
     geometry/Point.hpp
     geometry/Ray.hpp
     geometry/Segment.hpp
     geometry/Triangle.hpp
     geometry/UniformGrid.hpp

     ## numerics
     numerics/Determinants.hpp
     numerics/NumericArray.hpp
     numerics/Vector.hpp

     ## operators
     operators/detail/intersection_impl.hpp

     operators/intersection.hpp
     operators/orientation.hpp
     operators/squared_distance.hpp
     operators/closest_point.hpp

     ## spatial acceleration
     spatial_acceleration/BVHTree.hpp
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
    DEPENDS_ON            common slic
    )
