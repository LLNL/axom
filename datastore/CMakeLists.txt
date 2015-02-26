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
    ${PROJECT_SOURCE_DIR}/datastore/DataStore.cpps
    ${PROJECT_SOURCE_DIR}/datastore/DatastoreInterface.cpp
    )


#
# make the library
#
make_library(datastore ${datastore_sources})
add_dependencies(datastore copy_headers)
install(TARGETS datastore DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the
# build tree. It depends on "${datastore_headers}" to ensure it will be
# executed when a new header file is added.
#
add_custom_target(copy_headers
 COMMAND ${CMAKE_COMMAND}
         -DHEADER_INCLUDES_DIRECTORY=${HEADER_INCLUDES_DIRECTORY}
         -DLIBHEADERS="${datastore_headers}"
         -P ${PROJECT_SOURCE_DIR}/cmake/copy_headers.cmake

 DEPENDS
    ${datastore_headers}

 WORKING_DIRECTORY
    ${PROJECT_SOURCE_DIR}

 COMMENT
    "copy headers"
 )