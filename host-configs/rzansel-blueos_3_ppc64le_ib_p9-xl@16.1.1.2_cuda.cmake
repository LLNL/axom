#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.21.1/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/umpire-2023.06.0-5pojmv2de7q5pimrh253hwamycx7jfl7;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/raja-2023.06.0-7et5l25lkvvjz2r772uabj7hrstg3y2z;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/camp-2023.06.0-q4hulo6jzt2wnuzzj5jjmlvkdtu4ckdj;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/cub-2.1.0-45sgqog6spfots46rohzcokfoaxjil4p;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/mfem-4.5.2-264o4xkxhcafenixttd6lw5k466a62b3;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/hypre-2.24.0-x2zforves4hmpcuwjskxbvxu55i5kpfx;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/lua-5.4.4-6yyyruucztcz5p4w33mhydqn5kiywkjn;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/conduit-0.8.8-pil6joltz2jsp4vcdsfgndntfkqoi3na;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/parmetis-4.0.3-ogq2nkkjvarhjgaweb6b5lglfcmes3hm;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/metis-5.1.0-7ojc3lx2y6dk2g7gvl6m5n4fq5yewukb;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/hdf5-1.8.22-y23h3hgemr7s6terq2mvsz5if5u327sr;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/zlib-1.2.13-nd6r3ccnozxmd7jdwyxbgvqg5y2bn3ds;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/c2c-1.8.0-ckwviu4jwtot52xykd3kuryozuw5o35p;/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2/blt-0.5.3-47zfyusaljp2e7x3lsq3iggqvapxox6t;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/python-3.10.10;/usr/tcetmp/packages/lapack/lapack-3.9.0-P9-xl-2020.11.12;/usr/tce/packages/clang/clang-10.0.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/graphviz-7.1.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/doxygen-1.9.6;/usr/tce/packages/cuda/cuda-11.2.0;/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/latest/cppcheck-2.9;/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2022.08.19-cuda-11.2.0;/usr/tcetmp;/usr/tce/packages/cmake/cmake-3.21.1" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: xl@=16.1.1.2
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/spack/lib/spack/env/xl/xlc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/spack/lib/spack/env/xl/xlc++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/spack/lib/spack/env/xl/xlf90" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/xl/xl-2022.08.19-cuda-11.2.0/bin/xlc" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/xl/xl-2022.08.19-cuda-11.2.0/bin/xlC" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/xl/xl-2022.08.19-cuda-11.2.0/bin/xlf2003" CACHE PATH "")

endif()

set(CMAKE_C_FLAGS "--gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1" CACHE STRING "")

set(CMAKE_CXX_FLAGS "--gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1" CACHE STRING "")

set(CMAKE_Fortran_FLAGS "--gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1" CACHE STRING "")

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2022.08.19-cuda-11.2.0/bin/mpixlc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2022.08.19-cuda-11.2.0/bin/mpixlC" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2022.08.19-cuda-11.2.0/bin/mpixlf" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-xl-2022.08.19-cuda-11.2.0/bin/mpirun" CACHE PATH "")

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

set(CMAKE_CUDA_FLAGS "-restrict --expt-extended-lambda -Xcompiler=--gcc-toolchain=/usr/tce/packages/gcc/gcc-8.3.1 " CACHE STRING "")

# nvcc does not like gtest's 'pthreads' flag

set(gtest_disable_pthreads ON CACHE BOOL "")

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP OFF CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS OFF CACHE BOOL "")

set(BLT_EXE_LINKER_FLAGS "${BLT_EXE_LINKER_FLAGS} -Wl,-rpath,/usr/tce/packages/xl/xl-2022.08.19-cuda-11.2.0/lib" CACHE STRING "Adds a missing rpath for libraries associated with the fortran compiler")

set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-rpath,/usr/tce/packages/xl/xl-2022.08.19-cuda-11.2.0/lib" CACHE STRING "Adds a missing rpath for libraries associated with the fortran compiler")

set(BLT_FORTRAN_FLAGS "-WF,-C!  -qxlf2003=polymorphic" CACHE STRING "Converts C-style comments to Fortran style in preprocessed files")

set(BLT_CMAKE_IMPLICIT_LINK_DIRECTORIES_EXCLUDE "/usr/tce/packages/gcc/gcc-4.9.3/lib64;/usr/tce/packages/gcc/gcc-4.9.3/lib64/gcc/powerpc64le-unknown-linux-gnu/4.9.3;/usr/tce/packages/gcc/gcc-4.9.3/gnu/lib64;/usr/tce/packages/gcc/gcc-4.9.3/gnu/lib64/gcc/powerpc64le-unknown-linux-gnu/4.9.3" CACHE STRING "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/blueos_3_ppc64le_ib_p9/2023_10_17_11_51_45/xl-16.1.1.2" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-pil6joltz2jsp4vcdsfgndntfkqoi3na" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-ckwviu4jwtot52xykd3kuryozuw5o35p" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-264o4xkxhcafenixttd6lw5k466a62b3" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-y23h3hgemr7s6terq2mvsz5if5u327sr" CACHE PATH "")

set(LUA_DIR "${TPL_ROOT}/lua-5.4.4-6yyyruucztcz5p4w33mhydqn5kiywkjn" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-7et5l25lkvvjz2r772uabj7hrstg3y2z" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-5pojmv2de7q5pimrh253hwamycx7jfl7" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-q4hulo6jzt2wnuzzj5jjmlvkdtu4ckdj" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/blueos_3_ppc64le_ib_p9/2023_10_17_10_41_04/._view/5gug4et2qymmckk7yvtub4qcapp3figm" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/usr/tce/packages/clang/clang-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "/collab/usr/gapps/shroud/public/blueos_3_ppc64le_ib_p9/shroud-0.13.0/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.9.6/bin/doxygen" CACHE PATH "")


