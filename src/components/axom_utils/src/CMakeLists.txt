#
# Specify all headers
#
set(axom_utils_headers
    utilities/FileUtilities.hpp
    utilities/StringUtilities.hpp
    utilities/Timer.hpp
    utilities/Utilities.hpp

    detail/ChronoTimer.hpp
    detail/TimeofdayTimer.hpp

    numerics/Determinants.hpp
    numerics/linear_solve.hpp
    numerics/LU.hpp
    numerics/Matrix.hpp
    )

#
# Specify all sources
#
set(axom_utils_sources
    utilities/FileUtilities.cpp
    utilities/StringUtilities.cpp
    utilities/Utilities.cpp
    )

#
# make the library
#
blt_add_library( NAME                  axom_utils
                 SOURCES              "${axom_utils_sources}"
                 HEADERS              "${axom_utils_headers}"
                 HEADERS_OUTPUT_SUBDIR axom_utils
                 )
