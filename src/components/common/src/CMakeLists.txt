#
# Specify all headers
#
set(common_headers
    ATKMacros.hpp
    CommonTypes.hpp
    Utilities.hpp
    )

#
# Specify all sources
#
set(common_sources
	Utilities.cpp
    )

#
# make the library
#
if (ENABLE_MPI)
    make_library( LIBRARY_NAME
                      common
                  LIBRARY_SOURCES
                      "${common_sources}"
                      "${common_headers}"
                  WITH_MPI )
else()
    make_library( LIBRARY_NAME
                      common
                  LIBRARY_SOURCES
                      "${common_sources}"
                      "${common_headers}" )
endif()

add_dependencies(common copy_headers_common)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( ${PROJECT_NAME} "${common_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )
