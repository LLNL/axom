# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

################################
# Check necessary dependencies
################################
axom_component_requires(NAME       PRIMAL
                        COMPONENTS SLIC)

################################
# Specify all headers/sources
################################
set( primal_headers

     ## geometry
     geometry/BoundingBox.hpp
     geometry/OrientedBoundingBox.hpp
     geometry/OrientationResult.hpp
     geometry/NumericArray.hpp
     geometry/Plane.hpp
     geometry/Point.hpp
     geometry/Polygon.hpp
     geometry/Ray.hpp
     geometry/Segment.hpp
     geometry/Sphere.hpp
     geometry/Tetrahedron.hpp
     geometry/Triangle.hpp
     geometry/Vector.hpp

     ## operators
     operators/clip.hpp
     operators/closest_point.hpp
     operators/intersect.hpp
     operators/orientation.hpp
     operators/squared_distance.hpp
     operators/compute_bounding_box.hpp
     operators/in_sphere.hpp

     operators/detail/clip_impl.hpp
     operators/detail/intersect_impl.hpp

     ## spatial acceleration
     spatial_acceleration/BVH.hpp

     spatial_acceleration/linear_bvh/BVHData.hpp
     spatial_acceleration/linear_bvh/RadixTree.hpp
     spatial_acceleration/linear_bvh/aabb.hpp
     spatial_acceleration/linear_bvh/emit_bvh.hpp
     spatial_acceleration/linear_bvh/build_radix_tree.hpp
     spatial_acceleration/linear_bvh/math.hpp
     spatial_acceleration/linear_bvh/morton_codes.hpp
     spatial_acceleration/linear_bvh/policies.hpp
     spatial_acceleration/linear_bvh/range.hpp
     spatial_acceleration/linear_bvh/vec.hpp
   )

set( primal_sources

     ## operators
     operators/detail/intersect_impl.cpp

     ## spatial acceleration
     spatial_acceleration/linear_bvh/morton_codes.cpp
   )

################################
# Make/Install the library
################################
set( primal_dependencies
     core
     slic
     )

blt_list_append( TO primal_dependencies ELEMENTS raja IF ${RAJA_FOUND} )
blt_list_append( TO primal_dependencies ELEMENTS umpire IF ${UMPIRE_FOUND} )
blt_list_append( TO primal_dependencies ELEMENTS openmp IF ${ENABLE_OPENMP} )

blt_add_library(
    NAME                  primal
    SOURCES               ${primal_sources}
    HEADERS               ${primal_headers}
    FOLDER                axom/primal
    DEPENDS_ON            ${primal_dependencies}
    OBJECT                TRUE
    )

axom_write_unified_header(NAME    primal
                          HEADERS ${primal_headers})

axom_install_component(NAME    primal
                       HEADERS ${primal_headers}
                       )

################################
# add examples
################################
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# add tests
################################
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# Add code checks
################################
axom_add_code_checks(PREFIX primal)
