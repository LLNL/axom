#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/packages/cmake/cmake-3.19.2/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/umpire-2022.03.1-iut74u3azw6skamlcn2vhiyldxaoc3ud;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/raja-2022.03.1-skzhrbaucvynsb52vj4cavlyrgurcgu4;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/camp-2022.10.1-bbgab7ltyl3ss2pxx4rur6yd3q4hipuf;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/mfem-4.5.2-do5ocbeuy2moa5cytcsucxnhhy6ree5m;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/hypre-2.24.0-bwh42bebp4hiuwzlwcthpgkawape6dp3;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/conduit-0.8.8-o2daxjhl4oygijft7yjsgsx4f2kq3ogj;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/parmetis-4.0.3-ucmsjhvhd2okxlbbgovbpoq7vigmycbt;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/metis-5.1.0-ctoim5b4oc5mghlvvjwnhc2hsyzq6dca;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/hdf5-1.8.22-ducec2x4ah3xwjs2nk4eilysbup65mic;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/c2c-1.8.0-ns3fs5uxezevm2b7rp6y4bxvltrlwzaa;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/gmake-4.4.1-6lgy76r7dbvqm6mbwb2jkpcuycf7tb27;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6/blt-0.5.3-a47mj62hqv66gff55zxsxwissyyy5zwm;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/llvm-10.0.0;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/doxygen-1.8.14;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/cppcheck-2.9;/usr/tce/packages/mvapich2/mvapich2-2.3.6-clang-14.0.6;/usr/tce/packages/cmake/cmake-3.19.2;/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6;/usr/tce/packages/clang/clang-14.0.6" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=14.0.6
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/spack/lib/spack/env/clang/gfortran" CACHE PATH "")

else()

  set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-14.0.6/bin/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-14.0.6/bin/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/tce/packages/gcc/gcc-10.3.1/bin/gfortran" CACHE PATH "")

endif()

set(CMAKE_GENERATOR "Unix Makefiles" CACHE STRING "")

set(ENABLE_FORTRAN ON CACHE BOOL "")

set(BLT_EXE_LINKER_FLAGS " -Wl,-rpath,/usr/tce/packages/clang/clang-14.0.6/lib" CACHE STRING "Adds a missing libstdc++ rpath")

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-clang-14.0.6/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-clang-14.0.6/bin/mpicxx" CACHE PATH "")

set(MPI_Fortran_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.6-clang-14.0.6/bin/mpif90" CACHE PATH "")

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

#------------------------------------------------------------------------------
# Hardware
#------------------------------------------------------------------------------

#------------------------------------------------
# Hardware Specifics
#------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "")

#------------------------------------------------------------------------------
# TPLs
#------------------------------------------------------------------------------

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_07_09_19_16_24/clang-14.0.6" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-o2daxjhl4oygijft7yjsgsx4f2kq3ogj" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-ns3fs5uxezevm2b7rp6y4bxvltrlwzaa" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-do5ocbeuy2moa5cytcsucxnhhy6ree5m" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-ducec2x4ah3xwjs2nk4eilysbup65mic" CACHE PATH "")

set(LUA_DIR "/usr" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2022.03.1-skzhrbaucvynsb52vj4cavlyrgurcgu4" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2022.03.1-iut74u3azw6skamlcn2vhiyldxaoc3ud" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2022.10.1-bbgab7ltyl3ss2pxx4rur6yd3q4hipuf" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/2023_05_18_11_52_05/._view/btoxy5ovdbouub2brzxcmjwzdhvzatlc" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/llvm-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.8.14/bin/doxygen" CACHE PATH "")


