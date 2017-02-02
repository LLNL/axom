
################################################################################
# This directory contains source code for utility executables that utilize
# the toolkit or that support the toolkit
################################################################################


##################################################
# convert_protocol is a utility that converts
# a datastore from the sidre_hdf5 protocol
# to any of the other supported protocols
##################################################
set(datastore_converter_sources
    convert_sidre_protocol.cpp
)
set(datastore_converter_depends 
    spio sidre conduit conduit_relay hdf5 
    slic fmt 
    slam 
    )

if(ENABLE_LUMBERJACK)
  list(APPEND datastore_converter_depends lumberjack)
endif()


if(ENABLE_BOOST)
  list(APPEND datastore_converter_depends boost)
endif()


if(ENABLE_SLAM AND ENABLE_SPIO )
    blt_add_executable(
        NAME convert_sidre_protocol
        SOURCES ${datastore_converter_sources}
        OUTPUT_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        DEPENDS_ON ${datastore_converter_depends} sparsehash
        )
endif()

