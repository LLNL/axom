# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Spin component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       SPIN
                        COMPONENTS PRIMAL SLIC SLAM)

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
# Specify Spin headers/sources that depend on RAJA and Umpire
#------------------------------------------------------------------------------
if ( RAJA_FOUND AND UMPIRE_FOUND )

  set( spin_headers
       ${spin_headers}
       BVH.hpp

       ## internal
       internal/linear_bvh/BVHData.hpp
       internal/linear_bvh/BVH_impl.hpp
       internal/linear_bvh/QueryAccessor.hpp
       internal/linear_bvh/RadixTree.hpp
       internal/linear_bvh/build_radix_tree.hpp
       internal/linear_bvh/bvh_traverse.hpp
       internal/linear_bvh/bvh_vtkio.hpp
       internal/linear_bvh/emit_bvh.hpp

      )

endif()

#------------------------------------------------------------------------------
# Specify spin dependencies
#------------------------------------------------------------------------------
set( spin_depends_on
     core slic primal slam )

blt_list_append( TO spin_depends_on ELEMENTS sparsehash IF ${SPARSEHASH_FOUND} )
blt_list_append( TO spin_depends_on ELEMENTS RAJA IF ${RAJA_FOUND} )
blt_list_append( TO spin_depends_on ELEMENTS umpire IF ${UMPIRE_FOUND} )
blt_list_append( TO spin_depends_on ELEMENTS openmp  IF ${ENABLE_OPENMP} )
blt_list_append( TO spin_depends_on ELEMENTS cuda IF ${ENABLE_CUDA} )

#------------------------------------------------------------------------------
# Make/Install the library
#------------------------------------------------------------------------------
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

#------------------------------------------------------------------------------
# add examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

#------------------------------------------------------------------------------
# add tests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

#------------------------------------------------------------------------------
# Add code checks
#------------------------------------------------------------------------------
axom_add_code_checks(PREFIX spin)
