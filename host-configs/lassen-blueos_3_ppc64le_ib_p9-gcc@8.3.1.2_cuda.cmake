#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: gcc@=8.3.1.2
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_05_19_09_25_23/spack/lib/spack/env/gcc/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_05_19_09_25_23/spack/lib/spack/env/gcc/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_05_19_09_25_23/spack/lib/spack/env/gcc/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gfortran" CACHE PATH "")

endif()

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1/bin/mpirun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-np" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

set(BLT_MPI_COMMAND_APPEND "mpibind" CACHE STRING "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# Cuda
#------------------------------------------------

set(CUDA_TOOLKIT_ROOT_DIR "/usr/tce/packages/cuda/cuda-11.2.0" CACHE PATH "")

set(CMAKE_CUDA_COMPILER "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc" CACHE PATH "")

set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}" CACHE PATH "")

set(ENABLE_CUDA ON CACHE BOOL "")

set(CMAKE_CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "")

set(AXOM_ENABLE_ANNOTATIONS ON CACHE BOOL "")

set(CMAKE_CUDA_ARCHITECTURES "70" CACHE STRING "")

set(CMAKE_CUDA_FLAGS "-restrict --expt-extended-lambda " CACHE STRING "")

# nvcc does not like gtest's 'pthreads' flag

set(gtest_disable_pthreads ON CACHE BOOL "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")

set(BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE "/usr/tce/packages/gcc/gcc-4.9.3/lib64;/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/powerpc64le-unknown-linux-gnu/4.9.3;/usr/tce/packages/gcc/gcc-4.9.3/gnu/lib64;/usr/tce/packages/gcc/gcc-4.9.3/gnu/lib64/gcc/powerpc64le-unknown-linux-gnu/4.9.3" CACHE STRING "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_05_19_09_25_23/gcc-8.3.1.2" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.6-ike4aerqmgzxsamyrcmayhanwguxtykj" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-di6qfdy7o6arbmmjwxnm6wbxtf7igr4m" CACHE PATH "")

# MFEM not built

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-lit3gjb5rnt7q34sfiaf6amtl26xgtvp" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-44af23yo5cwnet2m4zcm2fsxhwfhjfdm" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.03.0-edfh2nhlthqqlki23gn7pbivjzi4ifmp" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.03.1-fz7kvvuvfgqhba6ehuwi5tndylwx4xbn" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.10.1-6dtococih64gzlco3d7g7skmaeqxsjv6" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/2023_04_18_13_41_48/._view/srxt35kojgk77f2222mk6mgv7z5jyyzz" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/usr/tce/packages/clang/clang-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.8.14/bin/doxygen" CACHE PATH "")


