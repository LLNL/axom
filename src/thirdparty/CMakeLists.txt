#
# Our third-party builtin libs are all currently header only,
# so here we are adding tests that make sure they link and can 
# be used
#

######################################################################
# Register the stla_io header-only library for reading ascii STL files
######################################################################
set(STLAIO_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/stla_io )
blt_register_library( 
    NAME stla_io
    INCLUDES ${STLAIO_INCLUDE_DIRS}
    )

############################################################
# Register the fmt header-only library for string formatting 
############################################################

# Note KW: Hack to copy fmt headers to proper directory.
#          Weird directory structure can be removed after using Chris White's blt fix for copy_headers macro
set(fmtHeaders thirdparty/fmt/fmt.hpp thirdparty/fmt/fmt/format.h thirdparty/fmt/fmt/format.cc)
set(fmtDest ${HEADER_INCLUDES_DIRECTORY}/fmt)
add_custom_target(copy_headers_fmt
    COMMAND             ${CMAKE_COMMAND} -DHEADER_INCLUDES_DIRECTORY=${fmtDest} -DLIBHEADERS="${fmtHeaders}" -P ${CMAKE_SOURCE_DIR}/blt/cmake/copy_headers.cmake
    WORKING_DIRECTORY   ${PROJECT_SOURCE_DIR}/..
    DEPENDS             custom_output
)
add_custom_command(OUTPUT custom_output COMMAND cmake -E echo)


blt_register_library( 
    NAME fmt
    INCLUDES ${HEADER_INCLUDES_DIRECTORY}
    )


################################
# Add tests
################################
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()


