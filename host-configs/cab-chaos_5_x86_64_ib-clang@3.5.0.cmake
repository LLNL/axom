##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# chaos_5_x86_64_ib-clang@3.5.0
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/cmake-3.3.1-yonvpqi2syehljvwq2uotcjmb2xydcoi/bin/cmake

#######
# using clang@3.5.0 compiler spec
#######

# c compiler used by spack
set("CMAKE_C_COMPILER" "/usr/global/tools/clang/chaos_5_x86_64_ib/clang-omp-3.5.0/bin/clang" CACHE PATH "")

# cpp compiler used by spack
set("CMAKE_CXX_COMPILER" "/usr/global/tools/clang/chaos_5_x86_64_ib/clang-omp-3.5.0/bin/clang++" CACHE PATH "")

# fortran compiler used by spack
# no fortran compiler

set("ENABLE_FORTRAN" "OFF" CACHE PATH "")

# hdf5 from uberenv
set("HDF5_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/hdf5-1.8.16-ffqxeq6flvcqk2fo6gscvkqvuhvouovx" CACHE PATH "")

# conduit from uberenv
set("CONDUIT_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/conduit-0.2.1-2zxz34weo5uxrn756dxq5dj2n2abzpvs" CACHE PATH "")

<<<<<<< HEAD:host-configs/cab-chaos_5_x86_64_ib-clang@3.5.0.cmake
# sparsehash headers from uberenv
set("SPARSEHASH_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/sparsehash-headers-2.0.2-afi2mzxg35xtfnrdfsci5toyjcwdolpq" CACHE PATH "")

=======
>>>>>>> 6fc779798a8ed375d0a6de251efd58f3ce4566f4:host-configs/surface-chaos_5_x86_64_ib-clang@3.5.0.cmake
# boost headers from uberenv
set("BOOST_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/boost-headers-1.58.0-a6n4lbbqdnhhkivpxazq5e4zrdhaspei" CACHE PATH "")

# python from uberenv
set("PYTHON_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/python-2.7.11-kmpji7fw4s22cfziy4byyfl22wjrmc7n/bin/python" CACHE PATH "")

# lua from uberenv
set("LUA_DIR" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/lua-5.1.5-mkz7ex36zg7oqdon6vxocr2tfk2i6in5" CACHE PATH "")

# doxygen from uberenv
set("DOXYGEN_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/doxygen-1.8.11-ef4b7mhn6rksthrqz47zxtndjhg2e6wq/bin/doxygen" CACHE PATH "")

# sphinx from uberenv
set("SPHINX_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/python-2.7.11-kmpji7fw4s22cfziy4byyfl22wjrmc7n/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set("UNCRUSTIFY_EXECUTABLE" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/uncrustify-0.61-hfvykjq7rsj5v7kb3blntogdcuewppl2/bin/uncrustify" CACHE PATH "")

# lcov and genhtml from uberenv
set("LCOV_PATH" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/lcov-1.11-otvm7uutci77rwj7cygg4d5vgecdorif/usr/bin/lcov" CACHE PATH "")

set("GENHTML_PATH" "/usr/workspace/wsa/axom/thirdparty_libs/builds/2017_03_22_09_30_19/spack/opt/spack/chaos_5_x86_64_ib/clang-3.5.0/lcov-1.11-otvm7uutci77rwj7cygg4d5vgecdorif/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc chaos5 gcc@4.7.1  host configs
##############################################################################

##############################################################################
# MPI - manually added for now.
##############################################################################
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/global/tools/clang/chaos_5_x86_64_ib/clang-omp-3.5.0/bin/mpiclang" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/global/tools/clang/chaos_5_x86_64_ib/clang-omp-3.5.0/bin/mpiclang++" CACHE PATH "")

set(MPIEXEC "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################

