# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

################################
# Check necessary dependencies
################################
axom_component_requires(NAME       SPIN
                        COMPONENTS PRIMAL SLIC SLAM)

################################
# Specify all headers/sources
################################
set( spin_headers

     BVHTree.hpp
     Brood.hpp
     DenseOctreeLevel.hpp
     ImplicitGrid.hpp
     MortonIndex.hpp
     OctreeBase.hpp
     OctreeLevel.hpp
     RectangularLattice.hpp
     SparseOctreeLevel.hpp
     SpatialOctree.hpp
     UniformGrid.hpp
   )

set( spin_sources
     ../Axom.cpp
   )

set( spin_depends_on
     core slic primal slam )

# Add sparsehash when available
if(SPARSEHASH_FOUND)
    list(APPEND spin_depends_on sparsehash)
endif()


################################
# Make/Install the library
################################
blt_add_library(
    NAME                  spin
    SOURCES               ${spin_sources}
    HEADERS               ${spin_headers}
    FOLDER                axom/spin
    DEPENDS_ON            ${spin_depends_on}
    OBJECT                TRUE )

axom_write_unified_header(NAME    spin
                          HEADERS ${spin_headers})

axom_install_component(NAME    spin
                       HEADERS ${spin_headers}
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
axom_add_code_checks(PREFIX spin)
