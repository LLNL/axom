
option(ENABLE_BOOST "Enables Boost" OFF)
option(ENABLE_PYTHON "Enables python use." ON)
option(ENABLE_CFORTRAN_API "Enables Fortran interface for components." ON)
option(ENABLE_ALL_COMPONENTS "Enables all components by default" ON)

if( (CMAKE_CXX_STANDARD EQUAL 11) OR (CMAKE_CXX_STANDARD EQUAL 14) )
	add_definitions("-DUSE_CXX11")
endif()
