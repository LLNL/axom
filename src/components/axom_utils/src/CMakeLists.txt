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

#
# Specify all headers
#
set(axom_utils_headers
    utilities/FileUtilities.hpp
    utilities/StringUtilities.hpp
    utilities/Timer.hpp
    utilities/Utilities.hpp
    utilities/About.hpp

    detail/ChronoTimer.hpp
    detail/TimeofdayTimer.hpp
    detail/TickCountTimer.hpp

    numerics/internal/matrix_norms.hpp

    numerics/Determinants.hpp
    numerics/LU.hpp
    numerics/Matrix.hpp
    numerics/eigen_solve.hpp
    numerics/eigen_sort.hpp
    numerics/linear_solve.hpp
    numerics/matvecops.hpp
    numerics/polynomial_solvers.hpp
    )

#
# Specify all sources
#
set(axom_utils_sources
    utilities/FileUtilities.cpp
    utilities/StringUtilities.cpp
    utilities/Utilities.cpp
    utilities/About.cpp

    numerics/polynomial_solvers.cpp
    )

#
# make the library
#
if (ENABLE_MPI)
    blt_add_library( NAME                  axom_utils
                     SOURCES              "${axom_utils_sources}"
                     HEADERS              "${axom_utils_headers}"
                     HEADERS_OUTPUT_SUBDIR axom_utils
                     DEPENDS_ON            mpi
                     FOLDER                axom/axom_utils
                     )
else()
    blt_add_library( NAME                  axom_utils
                     SOURCES              "${axom_utils_sources}"
                     HEADERS              "${axom_utils_headers}"
                     HEADERS_OUTPUT_SUBDIR axom_utils
                     FOLDER                axom/axom_utils
                     )

endif()
install(TARGETS              axom_utils
        EXPORT               axom_utils-targets
        DESTINATION          lib
        INCLUDES DESTINATION include/axom_utils
        )
install(FILES       ${axom_utils_headers}
        DESTINATION include/axom_utils
        )
