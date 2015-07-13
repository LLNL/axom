#######
# uberenv host-config for asctoolkit
#######
# cmake from uberenv
# cmake exectuable path: /home/taylor16/tpl/v1/spack/opt/x86_64/gcc@4.9.0/cmake@3.2.2/bin/cmake

# python from uberenv
set(PYTHON_EXECUTABLE "/home/taylor16/tpl/v1/spack/opt/x86_64/gcc@4.9.0/python@2.7.8/bin/python" CACHE PATH "")

# sphinx from uberenv
set(SPHINX_EXECUTABLE "/home/taylor16/tpl/v1/spack/opt/x86_64/gcc@4.9.0/python@2.7.8/bin/sphinx-build" CACHE PATH "")

# uncrustify from uberenv
set(UNCRUSTIFY_EXECUTABLE "/home/taylor16/tpl/v1/spack/opt/x86_64/gcc@4.9.0/uncrustify@0.61/bin/uncrustify" CACHE PATH "")

# boost headers from uberenv
set(ENABLE_BOOST ON CACHE PATH "")
set(BOOST_ROOT "/home/taylor16/tpl/v1/spack/opt/x86_64/gcc@4.9.0/boost-headers@1.58.0" CACHE PATH "")
#set(BOOST_ROOT "/home/taylor16/local/boost-1.58" CACHE PATH "")
