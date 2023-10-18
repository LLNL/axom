#------------------------------------------------------------------------------
# !!!! This is a generated file, edit at own risk !!!!
#------------------------------------------------------------------------------
# CMake executable path: /usr/tce/bin/cmake
#------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/umpire-2023.06.0-ubadmsj3l2zavunhkkt7nxevkpmssdvy;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/raja-2023.06.0-55ypap77cpmjgq24vaquanhsshjm3gqh;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/camp-2023.06.0-e6mknahsrmbbifbfv3umemetvqzoh3he;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/mfem-4.5.2-lwet4c3edrdvipij3r4itukmse2jklvy;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/hypre-2.24.0-bwh42bebp4hiuwzlwcthpgkawape6dp3;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/conduit-0.8.8-inqrhwboacvngevqbvavhpisx4aeqpvy;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/parmetis-4.0.3-jthqwflsj5bpe6onwsb5r2zi5wbx6pq4;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/metis-5.1.0-imzig3eih3itpztvvk4yildzgq6aeelo;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/hdf5-1.8.22-5mfm2r7fo5krpj7vvmayz24qksnluryl;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/c2c-1.8.0-5smd5ktj7yjc2egeyum4t5vlwmw3jktt;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/gmake-4.4.1-6lgy76r7dbvqm6mbwb2jkpcuycf7tb27;/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6/blt-0.5.3-xa3bmu75nkolqiewb5ftk5tp2jpzw57n;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/python-3.10.10;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/llvm-10.0.0;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/doxygen-1.9.6;/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/cppcheck-2.9;/usr/tce/packages/mvapich2/mvapich2-2.3.6-clang-14.0.6;/usr/tce;/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6;/usr/tce/packages/clang/clang-14.0.6" CACHE PATH "")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "")

#------------------------------------------------------------------------------
# Compilers
#------------------------------------------------------------------------------
# Compiler Spec: clang@=14.0.6
#------------------------------------------------------------------------------
if(DEFINED ENV{SPACK_CC})

  set(CMAKE_C_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/spack/lib/spack/env/clang/clang" CACHE PATH "")

  set(CMAKE_CXX_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/spack/lib/spack/env/clang/clang++" CACHE PATH "")

  set(CMAKE_Fortran_COMPILER "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/spack/lib/spack/env/clang/gfortran" CACHE PATH "")

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

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

set(ENABLE_MPI ON CACHE BOOL "")

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")

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

set(TPL_ROOT "/usr/WS1/axom/libs/toss_4_x86_64_ib/2023_10_18_08_27_25/clang-14.0.6" CACHE PATH "")

set(CONDUIT_DIR "${TPL_ROOT}/conduit-0.8.8-inqrhwboacvngevqbvavhpisx4aeqpvy" CACHE PATH "")

set(C2C_DIR "${TPL_ROOT}/c2c-1.8.0-5smd5ktj7yjc2egeyum4t5vlwmw3jktt" CACHE PATH "")

set(MFEM_DIR "${TPL_ROOT}/mfem-4.5.2-lwet4c3edrdvipij3r4itukmse2jklvy" CACHE PATH "")

set(HDF5_DIR "${TPL_ROOT}/hdf5-1.8.22-5mfm2r7fo5krpj7vvmayz24qksnluryl" CACHE PATH "")

set(LUA_DIR "/usr" CACHE PATH "")

set(RAJA_DIR "${TPL_ROOT}/raja-2023.06.0-55ypap77cpmjgq24vaquanhsshjm3gqh" CACHE PATH "")

set(UMPIRE_DIR "${TPL_ROOT}/umpire-2023.06.0-ubadmsj3l2zavunhkkt7nxevkpmssdvy" CACHE PATH "")

set(CAMP_DIR "${TPL_ROOT}/camp-2023.06.0-e6mknahsrmbbifbfv3umemetvqzoh3he" CACHE PATH "")

# scr not built

#------------------------------------------------------------------------------
# Devtools
#------------------------------------------------------------------------------

set(DEVTOOLS_ROOT "/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/2023_10_17_16_15_32/._view/ypfidjpdm7fhkqcpekza67w5xgiaw7yg" CACHE PATH "")

set(CLANGFORMAT_EXECUTABLE "/collab/usr/gapps/axom/devtools/toss_4_x86_64_ib/latest/llvm-10.0.0/bin/clang-format" CACHE PATH "")

set(PYTHON_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/python3.10" CACHE PATH "")

set(ENABLE_DOCS ON CACHE BOOL "")

set(SPHINX_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/sphinx-build" CACHE PATH "")

set(SHROUD_EXECUTABLE "${DEVTOOLS_ROOT}/python-3.10.10/bin/shroud" CACHE PATH "")

set(CPPCHECK_EXECUTABLE "${DEVTOOLS_ROOT}/cppcheck-2.9/bin/cppcheck" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "${DEVTOOLS_ROOT}/doxygen-1.9.6/bin/doxygen" CACHE PATH "")


