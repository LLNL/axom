#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# SYS_TYPE: blueos_3_ppc64le_ib_p9
# Compiler Spec: clang@upstream_xlf
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.14.5/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-upstream-2019.08.15/bin/clang" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-upstream-2019.08.15/bin/clang++" CACHE PATH "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER "/usr/tce/packages/xl/xl-2019.12.23/bin/xlf2003" CACHE PATH "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

# Root directory for generated TPLs
set(TPL_ROOT "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2020_06_16_18_16_44/clang-upstream_xlf" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.5.1" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.1.0" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.21" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.3.5" CACHE PATH "")

# SCR not built

set(RAJA_DIR "${TPL_ROOT}/raja-0.11.0" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2.1.0" CACHE PATH "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-upstream-2019.08.15/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-upstream-2019.08.15/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-upstream-2019.08.15/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-upstream-2019.08.15/bin/mpirun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-np" CACHE PATH "")

set(BLT_MPI_COMMAND_APPEND "mpibind" CACHE PATH "")

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# Root directory for generated developer tools
set(DEVTOOLS_ROOT "/usr/WS1/axom/devtools/blueos_3_ppc64le_ib_p9/2020_05_04_17_40_12/gcc-8.3.1" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.7/bin/python3.7" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.8.14/bin/doxygen" CACHE PATH "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.7/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.7/bin/shroud" CACHE PATH "")

set(UNCRUSTIFY_EXECUTABLE "${DEVTOOLS_ROOT}/uncrustify-0.61/bin/uncrustify" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-1.87/bin/cppcheck" CACHE PATH "")

#------------------------------------------------------------------------------
# Other machine specifics
#------------------------------------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

set(CMAKE_Fortran_COMPILER_ID "XL" CACHE PATH "Override to proper compiler family for XL")

set(BLT_FORTRAN_FLAGS "-WF,-C!  -qxlf2003=polymorphic" CACHE PATH "Converts C-style comments to Fortran style in preprocessed files")

set(BLT_EXE_LINKER_FLAGS "${BLT_EXE_LINKER_FLAGS} -Wl,-rpath,/usr/tce/packages/xl/xl-2019.12.23/lib" CACHE PATH "Adds a missing rpath for libraries associated with the fortran compiler")


