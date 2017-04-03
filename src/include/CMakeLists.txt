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
set(blt_dir ${CMAKE_SOURCE_DIR}/blt/cmake)

## custom target to copy the include headers
add_custom_target( copy_axom_includes ALL
        COMMAND ${CMAKE_COMMAND}
                 -DOUTPUT_DIRECTORY=${build_axom_includes}
                 -DLIBHEADERS="${axom_headers}"
                 -P ${blt_dir}/copy_headers.cmake

        DEPENDS
            ${arg_HEADERS}

        WORKING_DIRECTORY
            ${CMAKE_SOURCE_DIR}/include

        COMMENT
            "copy axom includes"
        )

include_directories( ${build_axom_includes} )

## add rules for the make install target
install( FILES ${axom_headers} DESTINATION ${install_axom_includes} )
