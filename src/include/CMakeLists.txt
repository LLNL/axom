###############################################################################
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and
# further review from Lawrence Livermore National Laboratory.
###############################################################################

## list of axom headers
set( axom_headers

     detail/cstdint_wrapper.hpp

     Macros.hpp
     Types.hpp
     )

## set include path for build and install spaces
set(build_axom_includes ${HEADER_INCLUDES_DIRECTORY}/axom )
set(install_axom_includes ${CMAKE_INSTALL_PREFIX}/include/axom )

## custom target to copy the include headers
add_custom_target( copy_axom_includes ALL
        COMMAND ${CMAKE_COMMAND}
                 -DOUTPUT_DIRECTORY=${build_axom_includes}
                 -DLIBHEADERS="${axom_headers}"
                 -P ${BLT_ROOT_DIR}/cmake/copy_headers.cmake

        DEPENDS
            ${axom_headers}

        WORKING_DIRECTORY
            ${CMAKE_SOURCE_DIR}/include

        COMMENT
            "copy axom includes"
        )

convert_to_native_escaped_file_path(${CMAKE_SOURCE_DIR} AXOM_SRC_DIR)
convert_to_native_escaped_file_path(${CMAKE_BINARY_DIR} AXOM_BIN_DIR)

configure_file(path_config.hpp.in "${build_axom_includes}/path_config.hpp")


include_directories( ${build_axom_includes} )

## add rules for the make install target
install( FILES ${axom_headers} DESTINATION ${install_axom_includes} )
install( FILES ${build_axom_includes}/config.hpp DESTINATION ${install_axom_includes} )
