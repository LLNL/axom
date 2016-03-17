#
# Specify all headers
#
set(common_headers
    ATKMacros.hpp
    CommonTypes.hpp
    FileUtilities.hpp
    StringUtilities.hpp
    Timer.hpp
    Utilities.hpp
    cstdint_wrapper.hpp
    
    detail/ChronoTimer.hpp
    detail/TimeofdayTimer.hpp
    )

#
# Specify all sources
#
set(common_sources
    FileUtilities.cpp
    StringUtilities.cpp
    Utilities.cpp
    )

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( ${PROJECT_NAME} "${common_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )

#
# make the library
#
make_library( LIBRARY_NAME
                  common
              LIBRARY_SOURCES
                  "${common_sources}"
                  "${common_headers}" )

install(TARGETS common DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

