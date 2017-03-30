#
# Specify all headers
#
set(axom_common_headers
    FileUtilities.hpp
    StringUtilities.hpp
    Timer.hpp
    Utilities.hpp

    detail/ChronoTimer.hpp
    detail/TimeofdayTimer.hpp
    )

#
# Specify all sources
#
set(axom_common_sources
    FileUtilities.cpp
    StringUtilities.cpp
    Utilities.cpp
    )

#
# make the library
#
blt_add_library( NAME                  axom_common
                 SOURCES              "${axom_common_sources}"
                 HEADERS              "${axom_common_headers}"
                 HEADERS_OUTPUT_SUBDIR axom_common
                 )
