#
# Specify all headers
#
set(common_headers
    Utilities.hpp
    Types.hpp
    )

#
# Specify all sources
#
set(common_sources
    )

#
# make the library
#
#make_library(common "${common_sources}" "${common_headers}")
#add_dependencies(common copy_headers_common)
#install(TARGETS common DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( ${PROJECT_NAME} "${common_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )
