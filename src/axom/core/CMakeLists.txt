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

################################
# Specify all headers/sources
################################
set(core_headers
    utilities/FileUtilities.hpp
    utilities/StringUtilities.hpp
    utilities/Timer.hpp
    utilities/Utilities.hpp
    utilities/About.hpp

    detail/cstdint_wrapper.hpp
    detail/ChronoTimer.hpp
    detail/TimeofdayTimer.hpp
    detail/TickCountTimer.hpp

    numerics/internal/matrix_norms.hpp

    numerics/Determinants.hpp
    numerics/jacobi_eigensolve.hpp
    numerics/LU.hpp
    numerics/Matrix.hpp
    numerics/eigen_solve.hpp
    numerics/eigen_sort.hpp
    numerics/linear_solve.hpp
    numerics/matvecops.hpp
    numerics/polynomial_solvers.hpp

    Array.hpp
    Macros.hpp
    Types.hpp
    )

set(core_sources
    utilities/FileUtilities.cpp
    utilities/StringUtilities.cpp
    utilities/Utilities.cpp
    utilities/About.cpp

    numerics/polynomial_solvers.cpp
    )

################################
# Make/Install the library
################################
if (ENABLE_MPI)
    blt_add_library( NAME                  core
                     SOURCES               ${core_sources}
                     HEADERS               ${core_headers}
                     DEPENDS_ON            mpi
                     FOLDER                axom/core
                     )
else()
    blt_add_library( NAME                  core
                     SOURCES               ${core_sources}
                     HEADERS               ${core_headers}
                     FOLDER                axom/core
                     )

endif()

axom_write_unified_header(NAME    core
                          HEADERS ${core_headers})

axom_install_component(NAME      core
                       HEADERS   ${core_headers}
                       )


################################
# Add tests
################################
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# Add docs
################################
if (AXOM_ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Add code checks
################################
axom_add_code_checks(PREFIX core)
