#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/umpire-2022.10.0-2hlqjzkv6skcbjnc2md2qj6n6ya52wny;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/raja-2022.10.5-ytnlbhcth3t4f5sbjgltkhldjb7ula5i;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/camp-2022.10.1-qhmylzgqn7rgm2dszeyxz4csc7ezztwb;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/cub-2.1.0-zpxcavwsek3eamis6c2ekyd4wpc2aff7;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/lua-5.4.4-uzsm5gerahtqcinqqfwleer62mlgl6xt;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/conduit-0.8.8-pslusbff27xjqw5sct3jnddyeqlyide3;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/parmetis-4.0.3-bcniwrmmygfl3ngqqy4ogyuopuilgfy2;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/metis-5.1.0-tht2h2hqk4psrvunj7egarksgyayew3d;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/hdf5-1.8.22-h2labsbvwxgxwvphjb5i357an3dgjhqp;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/zlib-1.2.13-obyznlavcvoyytpvwwyi65s5ngcnpigv;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/c2c-1.8.0-snpbqgmnulnruetig5di7nef2gir354f;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2/blt-0.5.3-umf5hhz3yl7vspds2fejsjdbmy2npyqb;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/usr/tce/packages/clang/clang-10.0.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/graphviz-7.1.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/doxygen-1.8.14;/usr/tce/packages/cuda/cuda-11.2.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/cppcheck-2.9;/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1;/usr/tcetmp;/usr/tce/packages/cmake/cmake-3.21.1" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: gcc@=8.3.1.2
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/spack/lib/spack/env/gcc/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/spack/lib/spack/env/gcc/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/spack/lib/spack/env/gcc/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gcc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/g++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gfortran" CACHE PATH "")

endif()

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

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

set(CUDAToolkit_ROOT "/usr/tce/packages/cuda/cuda-11.2.0" CACHE PATH "")

set(CMAKE_CUDA_COMPILER "${CUDAToolkit_ROOT}/bin/nvcc" CACHE PATH "")

set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}" CACHE PATH "")

set(CUDA_TOOLKIT_ROOT_DIR "/usr/tce/packages/cuda/cuda-11.2.0" CACHE PATH "")

set(CMAKE_CUDA_ARCHITECTURES "70" CACHE STRING "")

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

set(TPL_ROOT "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_07_09_16_57_46/gcc-8.3.1.2" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-pslusbff27xjqw5sct3jnddyeqlyide3" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-snpbqgmnulnruetig5di7nef2gir354f" CACHE PATH "")

# MFEM not built

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-h2labsbvwxgxwvphjb5i357an3dgjhqp" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-uzsm5gerahtqcinqqfwleer62mlgl6xt" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.10.5-ytnlbhcth3t4f5sbjgltkhldjb7ula5i" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.10.0-2hlqjzkv6skcbjnc2md2qj6n6ya52wny" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.10.1-qhmylzgqn7rgm2dszeyxz4csc7ezztwb" CACHE PATH "")

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


