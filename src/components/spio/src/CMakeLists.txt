#
# Specify all headers
#
set(spio_headers
    IOBaton.hpp
    IOManager.hpp
    )

#
# Specify all sources
#
set(spio_sources
    IOBaton.cpp
    IOManager.cpp
    )


#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target(  ${PROJECT_NAME}  "${spio_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )

#
# make the library
#
make_library( LIBRARY_NAME
                spio
              LIBRARY_SOURCES
                "${spio_sources}"
                "${spio_headers}"
              DEPENDS_ON
                conduit
                sidre
              WITH_MPI
            )

install(TARGETS spio DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
