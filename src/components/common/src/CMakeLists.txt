#
# Specify all headers
#
set(common_headers
    AxomMacros.hpp
    CommonTypes.hpp
    FileUtilities.hpp
    StringUtilities.hpp
    Timer.hpp
    Utilities.hpp
    
    detail/ChronoTimer.hpp
    detail/TimeofdayTimer.hpp
    detail/cstdint_wrapper.hpp
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
# make the library
#
blt_add_library( NAME
                    common
                 SOURCES
                    "${common_sources}"
                 HEADERS
                    "${common_headers}"
                 HEADERS_OUTPUT_SUBDIR
                    common )
