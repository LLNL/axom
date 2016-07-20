#
# CMakeConfigueFile.cmake - Create header of configuration options
#

if( (CMAKE_CXX_STANDARD EQUAL 11) OR (CMAKE_CXX_STANDARD EQUAL 14) )
    set(USE_CXX11 TRUE)
endif()

configure_file(
    components/common/src/config.hpp.in
    ${HEADER_INCLUDES_DIRECTORY}/common/config.hpp
)
