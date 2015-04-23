#
# Specify all headers
#
set(datastore_headers
    DataBuffer.hpp
    DataGroup.hpp
    DataStore.hpp
    DataView.hpp
    Utilities.hpp
    Types.hpp
    datastore.h
    )

#
# Specify all sources
#
set(datastore_sources
    DataBuffer.cpp
    DataGroup.cpp
    DataStore.cpp
    DataView.cpp
    Utilities.cpp
    )


#
# make the library
#
make_library(datastore "${datastore_sources}" "${datastore_headers}")
add_dependencies(datastore copy_headers_datastore)
add_dependencies(datastore conduit)
install(TARGETS datastore DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( ${PROJECT_NAME} "${datastore_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )
