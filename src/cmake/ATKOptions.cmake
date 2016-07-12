
option(ENABLE_BOOST "Enables Boost" OFF)
option(ENABLE_PYTHON "Enables python use." ON)
option(ENABLE_CFORTRAN_API "Enables Fortran interface for components." ON)
option(ENABLE_ALL_COMPONENTS "Enables all components by default" ON)
option(ENABLE_CXX11 "Enables C++11 if the compiler supports it" ON)

if ( NOT ENABLE_CXX11 )
    # do not use CXX11
elseif( NOT ( CMAKE_CXX_STANDARD LESS 11) )
    add_definitions("-DUSE_CXX11")
endif()