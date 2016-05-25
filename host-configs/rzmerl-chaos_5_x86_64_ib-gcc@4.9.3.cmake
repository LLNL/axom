##################################
# uberenv host-config
#
# This is a generated file, edit at own risk.
##################################
# chaos_5_x86_64_ib-gcc@4.9.3
##################################

# cmake from uberenv
# cmake executable path: /usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/cmake-3.3.1-3gc4unffj5rqcq35gecg4wv3roecpldt/bin/cmake

#######
# using gcc@4.9.3 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/apps/gnu/4.9.3/bin/gcc" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/apps/gnu/4.9.3/bin/g++" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/usr/apps/gnu/4.9.3/bin/gfortran" CACHE PATH "")

# hdf5 from uberenv
set(HDF5_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/hdf5-1.8.16-d3i3gs3i5iph3kakmjwen3eaiexwdrdg" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/conduit-github-2016-05-18-xqbkgfstnxnbt43ptpb6d26iv5pvytyk" CACHE PATH "")

# doxygen from uberenv
set(DOXYGEN_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/doxygen-1.8.11-7uy6qfezgsdt4dwwgrv6lkcn4g57brwo/bin/doxygen" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/python-2.7.11-nsotk32zeiza3u6gpkslyalrgfjmpspp/bin/python" CACHE PATH "")

# lua from uberenv
set(LUA_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/lua-5.1.5-axspkephlc7v6bqv5jnxcx4ronj57cd4" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/python-2.7.11-nsotk32zeiza3u6gpkslyalrgfjmpspp/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/uncrustify-0.61-px2meiscmkbwcnmmom3qnlzdzmf2yx7x/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/sparsehash-headers-2.0.2-w2nuzsfo46tzdm6kvghuzs4gaj7llnva" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/boost-headers-1.58.0-6wpeycxvyxiyu6g47hgwjefq4x5yur7w" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/lcov-1.11-dc44m4x5flevrfvqwytiscdc3dody4zx/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/2016_05_24_11_13_47/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/lcov-1.11-dc44m4x5flevrfvqwytiscdc3dody4zx/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
# Options added manually to 
# lc chaos5 gcc@4.9.3  host configs
##############################################################################

##############################################################################
# MPI - manually added for now
##############################################################################
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/local/tools/mvapich2-gnu-2.0/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/tools/mvapich2-gnu-2.0/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "/usr/local/tools/mvapich2-gnu-2.0/bin/mpif90" CACHE PATH "")

##############################################################################
# GCOV - manually added for now
##############################################################################
set(GCOV_PATH "/usr/apps/gnu/4.9.3/bin/gcov" CACHE PATH "")

##############################################################################
# !---------------------------------------------------------------------------
##############################################################################
