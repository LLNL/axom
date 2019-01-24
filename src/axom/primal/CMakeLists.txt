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
     geometry/MortonIndex.hpp
     geometry/NumericArray.hpp
     geometry/Plane.hpp
     geometry/Point.hpp
     geometry/Polygon.hpp
     geometry/Ray.hpp
     geometry/RectangularLattice.hpp
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
     spatial_acceleration/BVHTree.hpp
     spatial_acceleration/UniformGrid.hpp
   )

set( primal_sources

     operators/detail/intersect_impl.cpp
   )

################################
# Make/Install the library
################################
blt_add_library(
    NAME                  primal
    SOURCES               ${primal_sources}
    HEADERS               ${primal_headers}
    FOLDER                axom/primal
    DEPENDS_ON            core slic
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
# add docs
################################
if (AXOM_ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Add code checks
################################
axom_add_code_checks(PREFIX primal)
