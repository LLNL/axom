# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Axom's Core component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify all headers/sources
#------------------------------------------------------------------------------
set(core_headers

    ## utilities
    utilities/AnnotationMacros.hpp
    utilities/FileUtilities.hpp
    utilities/StringUtilities.hpp
    utilities/Timer.hpp
    utilities/Utilities.hpp
    utilities/About.hpp

    utilities/nvtx/interface.hpp
    utilities/nvtx/Macros.hpp
    utilities/nvtx/Range.hpp

    ## numerics
    numerics/internal/matrix_norms.hpp

    numerics/Determinants.hpp
    numerics/LU.hpp
    numerics/Matrix.hpp
    numerics/eigen_solve.hpp
    numerics/eigen_sort.hpp
    numerics/floating_point_limits.hpp
    numerics/jacobi_eigensolve.hpp
    numerics/linear_solve.hpp
    numerics/matvecops.hpp
    numerics/polynomial_solvers.hpp

    ## core
    Array.hpp
    Macros.hpp
    Path.hpp
    StackArray.hpp
    Types.hpp
    memory_management.hpp

    ## execution
    execution/execution_space.hpp
    execution/for_all.hpp
    execution/synchronize.hpp

    execution/internal/seq_exec.hpp
    execution/internal/omp_exec.hpp
    execution/internal/cuda_exec.hpp

    )

set(core_sources
    utilities/FileUtilities.cpp
    utilities/StringUtilities.cpp
    utilities/Utilities.cpp
    utilities/About.cpp

    utilities/nvtx/interface.cpp
    utilities/nvtx/Range.cpp

    numerics/polynomial_solvers.cpp

    Path.cpp
    Types.cpp
    )

#------------------------------------------------------------------------------
# Set library dependencies
#------------------------------------------------------------------------------

set( core_depends fmt )
blt_list_append( TO core_depends ELEMENTS umpire IF ${UMPIRE_FOUND} )
blt_list_append( TO core_depends ELEMENTS RAJA IF ${RAJA_FOUND} )
blt_list_append( TO core_depends ELEMENTS cuda nvToolsExt IF ${ENABLE_CUDA} )
blt_list_append( TO core_depends ELEMENTS openmp IF ${ENABLE_OPENMP} )
blt_list_append( TO core_depends ELEMENTS mpi IF ${ENABLE_MPI} )

#------------------------------------------------------------------------------
# Make/Install the library
#------------------------------------------------------------------------------

blt_add_library( NAME        core
                 SOURCES     ${core_sources}
                 HEADERS     ${core_headers}
                 DEPENDS_ON  ${core_depends}
                 FOLDER      axom/core
                 OBJECT      TRUE
                 )

axom_write_unified_header( NAME    core
                           HEADERS ${core_headers}
                           )

axom_install_component( NAME      core
                        HEADERS   ${core_headers}
                        )


#------------------------------------------------------------------------------
# Add tests
#------------------------------------------------------------------------------
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if(AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

#------------------------------------------------------------------------------
# Add code checks
#------------------------------------------------------------------------------
axom_add_code_checks(PREFIX core)
