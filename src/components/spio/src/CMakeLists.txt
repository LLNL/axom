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
# make the library
#
blt_add_library( NAME
                     spio
                 SOURCES
                     "${spio_sources}"
                 HEADERS
                     "${spio_headers}"
                 DEPENDS_ON
                     conduit
                     hdf5
                     sidre
                )

install(TARGETS spio DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
