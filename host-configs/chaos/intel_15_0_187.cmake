###############################################################################
#
# CMake Cache Seed file for chaos_5_x86_64_ib machines using icc 15 compiler
#
###############################################################################

##################################
# uberenv host-config
##################################
# chaos_5_x86_64_ib-intel@15.0.0
##################################

# cmake from uberenv
# cmake exectuable path: /usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.0/cmake-3.2.2-qjlmocqnctzvaswr6r5ftesbtofsp7yy/bin/cmake

#######
# using intel@15.0.0 compiler spec
#######

# c compiler used by spack
set(CMAKE_C_COMPILER "/usr/local/bin/icc-15.0.090" CACHE PATH "")

# cpp compiler used by spack
set(CMAKE_CXX_COMPILER "/usr/local/bin/icpc-15.0.090" CACHE PATH "")

# fortran compiler used by spack
set(ENABLE_FORTRAN ON CACHE PATH "")

set(CMAKE_Fortran_COMPILER  "/usr/local/bin/ifort-15.0.090" CACHE PATH "")

# conduit from uberenv
set(CONDUIT_DIR "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.0/conduit-github-ijpvf7vcibvjjszq2d34fe4ccrohdiqh" CACHE PATH "")

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.0/python-2.7.8-s2musnt6jmsyrcj2dsddvzirshcpjcbz/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.0/python-2.7.8-s2musnt6jmsyrcj2dsddvzirshcpjcbz/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.0/uncrustify-0.61-bjkdysaosc3zdzkveq2r4offg5iaa26h/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.0/sparsehash-headers-2.0.2-xlswjzupyehycxjk2rwrwb325ykw5bux" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.0/boost-headers-1.58.0-vr4eozxiybr54uz54ovcx2ghwueqjov4" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.0/lcov-1.11-fvr3ipjnzk4mnfghkhdblakjjlhlk7ep/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/spack/chaos_5_x86_64_ib/intel-15.0.0/lcov-1.11-fvr3ipjnzk4mnfghkhdblakjjlhlk7ep/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################


###############################################################################
# Additional Compiler Flags
###############################################################################

# additional flags for all configurations 
#set(EXTRA_CXX_FLAGS "-DEXTRA_FLAGS_CXX_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_DEBUG_FLAGS "-DEXTRA_CXX_DEBUG_FLAGS_DEFINE" CACHE PATH "")

# additional flags for debug builds
# set(EXTRA_CXX_RELEASE_FLAGS "-DEXTRA_CXX_RELEASE_FLAGS_DEFINE" CACHE PATH "")
