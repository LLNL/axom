# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Primal -- Computational geometry primitives
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       PRIMAL
                        COMPONENTS SLIC)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set( primal_headers

     ## geometry
     geometry/BezierCurve.hpp
     geometry/BoundingBox.hpp
     geometry/CurvedPolygon.hpp
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
     operators/detail/intersect_bezier_impl.hpp
     operators/detail/intersect_ray_impl.hpp
     operators/detail/intersect_bounding_box_impl.hpp
     operators/detail/intersect_impl.hpp
   )

set( primal_sources

     ## operators
     operators/detail/intersect_impl.cpp
   )

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set( primal_dependencies
     core
     slic
     )

blt_list_append( TO primal_dependencies ELEMENTS cuda IF ENABLE_CUDA )

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

#------------------------------------------------------------------------------
# Add tests and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

#------------------------------------------------------------------------------
# Add code checks
#------------------------------------------------------------------------------
axom_add_code_checks(PREFIX primal)
