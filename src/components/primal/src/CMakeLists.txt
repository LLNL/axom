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
     geometry/OrientedBoundingBox.hpp
     geometry/HyperSphere.hpp
     geometry/MortonIndex.hpp
     geometry/NumericArray.hpp
     geometry/Point.hpp
     geometry/Polygon.hpp
     geometry/Ray.hpp
     geometry/RectangularLattice.hpp
     geometry/Segment.hpp
     geometry/Triangle.hpp
     geometry/Vector.hpp

     ## operators
     operators/clip.hpp
     operators/closest_point.hpp
     operators/intersect.hpp
     operators/orientation.hpp
     operators/squared_distance.hpp
     operators/compute_bounding_box.hpp

     operators/detail/clip_impl.hpp
     operators/detail/intersect_impl.hpp

     ## spatial acceleration
     spatial_acceleration/BVHTree.hpp
     spatial_acceleration/UniformGrid.hpp
   )

# For now header only libraries cannot have dependencies until blt_setup_target is fixed
#   to toggle to INTERFACE only for them
#     DEPENDS_ON            slic
blt_add_library(
    NAME                  primal
    HEADERS               ${primal_headers}
    HEADERS_OUTPUT_SUBDIR primal
    )
install(TARGETS              primal
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include/primal
        )
install(FILES       ${primal_headers}
        DESTINATION include/primal
        )
