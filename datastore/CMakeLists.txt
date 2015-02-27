#
# Specify all headers
#
set(datastore_headers
    ${PROJECT_SOURCE_DIR}/datastore/Attributes.hpp
    ${PROJECT_SOURCE_DIR}/datastore/DataGroup.hpp
    ${PROJECT_SOURCE_DIR}/datastore/DataObject.hpp
    ${PROJECT_SOURCE_DIR}/datastore/DataStore.hpp
    ${PROJECT_SOURCE_DIR}/datastore/DatastoreInterface.hpp
    ${PROJECT_SOURCE_DIR}/datastore/Types.hpp
    ${PROJECT_SOURCE_DIR}/datastore/Utilities.hpp
    ${PROJECT_SOURCE_DIR}/datastore/datastore.h
    )

#
# Specify all sourcess
#
set(datastore_sources
    ${PROJECT_SOURCE_DIR}/datastore/DataGroup.cpp
    ${PROJECT_SOURCE_DIR}/datastore/DataObject.cpp
    ${PROJECT_SOURCE_DIR}/datastore/DataStore.cpp
    ${PROJECT_SOURCE_DIR}/datastore/DatastoreInterface.cpp
    )


#
# make the library
#
make_library(datastore "${datastore_sources}")
add_dependencies(datastore copy_headers)
install(TARGETS datastore DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( "${datastore_headers}" ${HEADER_INCLUDES_DIRECTORY} )