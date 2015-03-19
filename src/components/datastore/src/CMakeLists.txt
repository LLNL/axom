#
# Specify all headers
#
set(datastore_headers
    DataBuffer.hpp
    DataGroup.hpp
    DataObject.hpp
    DataStore.hpp
    DataView.hpp
    Types.hpp
    Utilities.hpp
    datastore.h
    )

#
# Specify all sourcess
#
set(datastore_sources
	DataBuffer.cpp
    DataGroup.cpp
    DataObject.cpp
    DataStore.cpp
	DataView.cpp
    )


#
# make the library
#
make_library(datastore "${datastore_sources}")
add_dependencies(datastore copy_headers)
add_dependencies(datastore conduit)
install(TARGETS datastore DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( "${datastore_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )