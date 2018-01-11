#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

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
            ${PROJECT_SOURCE_DIR}/src/include

        COMMENT
            "copy axom includes"
        )


include_directories( ${build_axom_includes} )

################################
# Add code checks
################################
axom_add_code_checks(PREFIX axom_includes)

## add rules for the make install target
install( FILES ${axom_headers} DESTINATION ${install_axom_includes} )
install( FILES ${build_axom_includes}/config.hpp DESTINATION ${install_axom_includes} )
