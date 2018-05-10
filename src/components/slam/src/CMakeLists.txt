#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

######################################
# Specify all headers and source files
######################################
set(slam_headers
    # Utility files
    Utilities.hpp
    FieldRegistry.hpp
    ModularInt.hpp

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
    #ProductSet.hpp
    DynamicSet.hpp
    Set.hpp
    BitSet.hpp
    internal/BitTwiddle.hpp

    # SRM Relation headers
    Relation.hpp
    StaticRelation.hpp
    DynamicVariableRelation.hpp

    # SRM Map headers
    Map.hpp
    )

set(slam_sources
    # SRM Set sources
    # RangeSet.cpp
    # IndirectionSet.cpp
    OrderedSet.cpp

    # SRM Relation sources
    Relation.cpp
    DynamicVariableRelation.cpp

    # SRM Map sources

    # Other sources
    Utilities.cpp
    )


#get_source_file_property(_origflags StaticVariableRelation.hpp COMPILE_FLAGS)
#set_source_files_properties(${slam_headers} "${_origflags} -fwhole-program -funroll-loops ")
#set_source_files_properties(${slam_sources} "${_origflags} -fwhole-program -funroll-loops ")


####################################
# make the library
####################################

set(slam_depends_on axom_utils slic)

blt_add_library(
    NAME        slam
    SOURCES     ${slam_sources}
    HEADERS     ${slam_headers}
    HEADERS_OUTPUT_SUBDIR slam
    DEPENDS_ON  ${slam_depends_on}
    FOLDER      axom/slam
    )
install(TARGETS              slam
        EXPORT               slam-targets
        DESTINATION          lib
        INCLUDES DESTINATION include/slam
        )
install(FILES       ${slam_headers}
        DESTINATION include/slam
        )

