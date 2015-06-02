#
# Specify all headers
#

add_subdirectory(capi)

set(sidre_headers
    core/Collections.hpp
    core/DataBuffer.hpp
    core/DataGroup.hpp
    core/DataStore.hpp
    core/DataView.hpp
    core/SidreTypes.hpp
    core/sidre.hpp
    ${capi_headers}
)

#
# Specify all sources
#
set(sidre_sources
    core/DataBuffer.cpp
    core/DataGroup.cpp
    core/DataStore.cpp
    core/DataView.cpp
    ${capi_source}
)

#
# make the library
#
make_library(sidre "${sidre_sources}")
add_dependencies(sidre sidre_generate)
add_dependencies(sidre copy_headers_sidre)
add_dependencies(sidre copy_headers_common)
add_dependencies(sidre copy_headers_conduit)
add_dependencies(sidre conduit)
install(TARGETS sidre DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( ${PROJECT_NAME} "${sidre_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )
add_dependencies(copy_headers_sidre sidre_generate)
