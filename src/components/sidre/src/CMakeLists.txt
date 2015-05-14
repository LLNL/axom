#
# Specify all headers
#
set(sidre_headers
    Collections.hpp
    DataBuffer.hpp
    DataGroup.hpp
    DataStore.hpp
    DataView.hpp
    sidre.hpp
    )

#
# Specify all sources
#
set(sidre_sources
    DataBuffer.cpp
    DataGroup.cpp
    DataStore.cpp
    DataView.cpp
    )


#
# make the library
#
make_library(sidre "${sidre_sources}" "${sidre_headers}")
add_dependencies(sidre copy_headers_sidre)
add_dependencies(sidre copy_headers_common)
add_dependencies(sidre conduit)
install(TARGETS sidre DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( ${PROJECT_NAME} "${sidre_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )
