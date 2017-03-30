#
# Specify all headers
#
set(axom_common_headers
    utilities/FileUtilities.hpp
    utilities/StringUtilities.hpp
    utilities/Timer.hpp
    utilities/Utilities.hpp

    detail/ChronoTimer.hpp
    detail/TimeofdayTimer.hpp
    )

#
# Specify all sources
#
set(axom_common_sources
    utilities/FileUtilities.cpp
    utilities/StringUtilities.cpp
    utilities/Utilities.cpp
    )

#
# make the library
#
blt_add_library( NAME                  axom_common
                 SOURCES              "${axom_common_sources}"
                 HEADERS              "${axom_common_headers}"
                 HEADERS_OUTPUT_SUBDIR axom_common
                 )
