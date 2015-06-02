#
# Specify all headers
#
set(meshapi_headers
    hello.hpp
    )

#
# Specify all sources
#
set(meshapi_sources
    hello.cpp
    )


#
# make the library
#
make_library( LIBRARY_NAME
                meshapi
              LIBRARY_SOURCES
                "${meshapi_sources}"
                "${meshapi_headers}"
            )

add_dependencies(meshapi copy_headers_meshapi)
add_dependencies(meshapi sidre)
install(TARGETS meshapi DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target(  ${PROJECT_NAME}  "${meshapi_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )
