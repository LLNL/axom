INCLUDE_DIRECTORIES(
    ${Boost_INCLUDE_DIRS}
)

######################################
# Specify all headers and source files
######################################
set(meshapi_headers
    # Utility files
    Utilities.hpp
    FileUtilities.hpp
    FieldRegistry.hpp
    
    #SRM policies
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
    Set.hpp
    
    # SRM Relation headers
    Relation.hpp
    StaticVariableRelation.hpp
    StaticConstantRelation.hpp
    DynamicVariableRelation.hpp
    
    # SRM Map headers
    Map.hpp
    )

set(meshapi_sources
    # SRM Set sources
    # RangeSet.cpp
    # IndirectionSet.cpp
    OrderedSet.cpp
    
    # SRM Relation sources
    Relation.cpp
    StaticVariableRelation.cpp
    #StaticConstantRelation.cpp
    DynamicVariableRelation.cpp
    
    # SRM Map sources
    )


#get_source_file_property(_origflags StaticVariableRelation.hpp COMPILE_FLAGS)
#set_source_files_properties(${meshapi_headers} "${_origflags} -fwhole-program -funroll-loops ")
#set_source_files_properties(${meshapi_sources} "${_origflags} -fwhole-program -funroll-loops ")


#########################################################################
# custom target to copy header files in the "include" folder of the build tree.
#########################################################################

copy_headers_target(  ${PROJECT_NAME}  "${meshapi_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )



####################################=
# make the library
####################################

make_library( LIBRARY_NAME
                meshapi
              LIBRARY_SOURCES
                "${meshapi_sources}"
                "${meshapi_headers}"
              DEPENDS_ON
                slic
                sidre
            )

####################################=
# install target for meshapi
####################################

install(TARGETS meshapi DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
