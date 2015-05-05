INCLUDE_DIRECTORIES(  
    ${Boost_INCLUDE_DIRS}   
)

#
# Specify all headers
#
set(meshapi_headers
    Hello.hpp
    Utilities.hpp
	# SRM Set headers
    Set.hpp
    OrderedSet.hpp
	# SRM Relation headers
	Relation.hpp
	StaticVariableRelation.hpp
	# SRM Map headers
    )

#
# Specify all sources
#
set(meshapi_sources
    Hello.cpp
    # SRM Set sources
    Set.cpp
    OrderedSet.cpp
    # SRM Relation sources
    StaticVariableRelation.cpp
    # SRM Map sources
    )


#
# make the library
#
make_library(meshapi "${meshapi_sources}" "${meshapi_headers}")
add_dependencies(meshapi copy_headers_meshapi)
add_dependencies(meshapi datastore)
install(TARGETS meshapi DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target(  ${PROJECT_NAME}  "${meshapi_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )
