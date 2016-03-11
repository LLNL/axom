##################################
# uberenv host-config
##################################
# chaos_5_x86_64_ib-gcc@4.9.3
##################################

# cmake from uberenv
# cmake executable path: /usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/cmake-3.2.2-v2wm3pqxrvtdqgl6ja6zixly3adz3m5y/bin/cmake

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

# conduit from uberenv
set(CONDUIT_DIR "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/conduit-github-gbnlaqxhiuu7sobrv4utggdnwwvl6llv" CACHE PATH "")

# doxygen from uberenv
# doxygen exectuable path: /usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/doxygen-1.8.10-auub5iinfei65n3qne2yiqkv24y2ww42/bin/doxygen

# python from uberenv
set(PYTHON_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/python-2.7.8-ed7poggszl4rmehk6cmd6d2gfqhfdrnn/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/python-2.7.8-ed7poggszl4rmehk6cmd6d2gfqhfdrnn/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/uncrustify-0.61-px2meiscmkbwcnmmom3qnlzdzmf2yx7x/bin/uncrustify" CACHE PATH "")

# sparsehash headers from uberenv
set(SPARSEHASH_DIR "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/sparsehash-headers-2.0.2-w2nuzsfo46tzdm6kvghuzs4gaj7llnva" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/boost-headers-1.58.0-6wpeycxvyxiyu6g47hgwjefq4x5yur7w" CACHE PATH "")

# lcov and genhtml from uberenv
set(LCOV_PATH "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/lcov-1.11-dc44m4x5flevrfvqwytiscdc3dody4zx/usr/bin/lcov" CACHE PATH "")

set(GENHTML_PATH "/usr/gapps/asctoolkit/thirdparty_libs/2016_03_11/spack/opt/spack/chaos_5_x86_64_ib/gcc-4.9.3/lcov-1.11-dc44m4x5flevrfvqwytiscdc3dody4zx/usr/bin/genhtml" CACHE PATH "")

##################################
# end uberenv host-config
##################################

