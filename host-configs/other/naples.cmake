###############################################################################
#
# CMake Cache Seed file for Cyrus' laptop.
#
###############################################################################

###############################################################################
# use clang compilers
###############################################################################
set(CMAKE_C_COMPILER "clang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "clang++" CACHE PATH "")

#######
# uberenv host-config for asctoolkit
#######
# cmake from uberenv
# cmake exectuable path: /Users/harrison37/Work/asctoolkit/uberenv_libs/spack/opt/macosx_10.9_x86_64/gcc@4.2.1/cmake@3.2.2/bin/cmake

# python from uberenv
set(PYTHON_EXECUTABLE "/Users/harrison37/Work/asctoolkit/uberenv_libs/spack/opt/macosx_10.9_x86_64/gcc@4.2.1/python@2.7.8/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/Users/harrison37/Work/asctoolkit/uberenv_libs/spack/opt/macosx_10.9_x86_64/gcc@4.2.1/python@2.7.8/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/Users/harrison37/Work/asctoolkit/uberenv_libs/spack/opt/macosx_10.9_x86_64/gcc@4.2.1/uncrustify@0.61/bin/uncrustify" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/Users/harrison37/Work/asctoolkit/uberenv_libs/spack/opt/macosx_10.9_x86_64/gcc@4.2.1/boost-headers@1.58.0" CACHE PATH "")



