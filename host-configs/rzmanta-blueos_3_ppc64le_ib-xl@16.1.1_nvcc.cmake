#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.14.5/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: xl@16.1.1_nvcc
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib/2021_06_08_16_11_03/spack/lib/spack/env/xl/xlc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib/2021_06_08_16_11_03/spack/lib/spack/env/xl/xlc++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib/2021_06_08_16_11_03/spack/lib/spack/env/xl/xlf90" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/xl/xl-2019.08.20/bin/xlc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/xl/xl-2019.08.20/bin/xlC" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/xl/xl-2019.08.20/bin/xlf2003" CACHE PATH "")

endif()

set(CMAKE_Fortran_COMPILER_ID "XL" CACHE STRING "Override to proper compiler family for XL")

set(CMAKE_C_COMPILER_ID "XL" CACHE STRING "Override to proper compiler family for XL")

set(CMAKE_CXX_COMPILER_ID "XL" CACHE STRING "Override to proper compiler family for XL")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2019.08.20/bin/mpixlc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2019.08.20/bin/mpixlC" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2019.08.20/bin/mpixlf" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2019.08.20/bin/mpirun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

set(BLT_MPI_COMMAND_APPEND "mpibind" CACHE STRING "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# Cuda
#------------------------------------------------

set(CUDA_TOOLKIT_ROOT_DIR "/usr/tce/packages/cuda/cuda-10.1.243" CACHE PATH "")

set(CMAKE_CUDA_COMPILER "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc" CACHE PATH "")

set(CMAKE_CUDA_HOST_COMPILER "${MPI_CXX_COMPILER}" CACHE PATH "")

set(ENABLE_CUDA ON CACHE BOOL "")

set(CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "")

set(AXOM_ENABLE_ANNOTATIONS ON CACHE BOOL "")

set(CMAKE_CUDA_ARCHITECTURES "60" CACHE STRING "")

set(CMAKE_CUDA_FLAGS "-restrict --expt-extended-lambda -arch sm_${CMAKE_CUDA_ARCHITECTURES}  -std=c++11" CACHE STRING "")

# nvcc does not like gtest's 'pthreads' flag

set(gtest_disable_pthreads ON CACHE BOOL "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")

set(BLT_FORTRAN_FLAGS "-WF,-C!  -qxlf2003=polymorphic" CACHE STRING "Converts C-style comments to Fortran style in preprocessed files")

set(BLT_EXE_LINKER_FLAGS "${BLT_EXE_LINKER_FLAGS} -Wl,-rpath,/usr/tce/packages/xl/xl-2019.08.20/lib" CACHE STRING "Adds a missing rpath for libraries associated with the fortran compiler")

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath,/usr/tce/packages/xl/xl-2019.08.20/lib" CACHE STRING "Adds a missing rpath for libraries associated with the fortran compiler")

set(BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE "/usr/tce/packages/gcc/gcc-4.9.3/lib64;/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/powerpc64le-unknown-linux-gnu/4.9.3" CACHE STRING "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

# Root directory for generated TPLs

set(TPL_ROOT "/usr/WS1/axom/libs/blueos_3_ppc64le_ib/2021_06_08_16_11_03/xl-16.1.1_nvcc" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.7.2" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.3.0" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.2.0" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.3.5" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-0.12.1" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-4.0.1" CACHE PATH "")

# scr not build

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

# Root directory for generated developer tools

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib/2020_08_21_20_29_39/gcc-8.3.1" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/usr/tce/packages/clang/clang-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.7/bin/python3.7" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.7/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.7.7/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-1.87/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.8.14/bin/doxygen" CACHE PATH "")


