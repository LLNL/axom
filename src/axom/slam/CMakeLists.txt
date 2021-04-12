# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# SLAM -- Set-theoretic Lightweight API for Meshes (low level mesh API)
#------------------------------------------------------------------------------
# Some additional flags that might be useful
#   To avoid warnings: -Wno-missing-braces
#   Other: -fno-strict-aliasing -Wextra
#   Vectorizing: -ftree-vectorize -fdump-tree-vect -ftree-vectorizer-verbose=7
#                -fdump-tree-optimized-details -march=native 
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Check necessary dependencies
#------------------------------------------------------------------------------
axom_component_requires(NAME       SLAM
                        COMPONENTS SLIC)


#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(slam_headers
    # Utility files
    Utilities.hpp
    FieldRegistry.hpp
    ModularInt.hpp
    IteratorBase.hpp

    #SRM policies
    policies/CardinalityPolicies.hpp
    policies/SizePolicies.hpp
    policies/OffsetPolicies.hpp
    policies/StridePolicies.hpp
    policies/IndirectionPolicies.hpp
    policies/SubsettingPolicies.hpp
    policies/PolicyTraits.hpp

    # SRM Set headers
    IndirectionSet.hpp
    NullSet.hpp
    OrderedSet.hpp
    RangeSet.hpp
    ProductSet.hpp
    DynamicSet.hpp
    Set.hpp
    BivariateSet.hpp
    BitSet.hpp
    internal/BitTwiddle.hpp
    RelationSet.hpp

    # SRM Relation headers
    Relation.hpp
    StaticRelation.hpp
    DynamicVariableRelation.hpp
    DynamicConstantRelation.hpp

    # SRM Map headers
    Map.hpp
    MapBase.hpp
    BivariateMap.hpp
    SubMap.hpp
    DynamicMap.hpp
    )

set(slam_sources
    # SRM Set sources
    # RangeSet.cpp
    # IndirectionSet.cpp
    BitSet.cpp
    OrderedSet.cpp

    # SRM Relation sources

    # SRM Map sources

    # Other sources
    )

#------------------------------------------------------------------------------
# Build and install the library
#------------------------------------------------------------------------------
set(slam_depends_on core slic)

blt_add_library(
    NAME        slam
    SOURCES     ${slam_sources}
    HEADERS     ${slam_headers}
    DEPENDS_ON  ${slam_depends_on}
    FOLDER      axom/slam
    OBJECT      TRUE )

axom_write_unified_header(NAME    slam
                          HEADERS ${slam_headers} )

axom_install_component(NAME    slam
                       HEADERS ${slam_headers} )

#------------------------------------------------------------------------------
# Add tests, benchmarks and examples
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
  if (ENABLE_BENCHMARKS)
    add_subdirectory(benchmarks)
  endif()
endif()

if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()


#------------------------------------------------------------------------------
# Add code checks
#------------------------------------------------------------------------------
axom_add_code_checks(
        PREFIX slam
        EXCLUDES
            slam/examples/lulesh2.0.3
            slam/examples/tinyHydro )

