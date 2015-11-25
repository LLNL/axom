#
# Specify all headers
#
set(ioparallel_headers
    IOBaton.hpp
    IOParallel.hpp
    )

#
# Specify all sources
#
set(ioparallel_sources
    IOBaton.cpp
    IOParallel.cpp
    )


#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target(  ${PROJECT_NAME}  "${ioparallel_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )

#
# make the library
#
make_library( LIBRARY_NAME
                ioparallel
              LIBRARY_SOURCES
                "${ioparallel_sources}"
                "${ioparallel_headers}"
              DEPENDS_ON
                conduit
                sidre
              WITH_MPI
            )

install(TARGETS ioparallel DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
