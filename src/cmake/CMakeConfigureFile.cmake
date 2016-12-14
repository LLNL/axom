#
# CMakeConfigueFile.cmake - Create header of configuration options
#

if( (CMAKE_CXX_STANDARD EQUAL 11) OR (CMAKE_CXX_STANDARD EQUAL 14) )
    set(ATK_USE_CXX11 TRUE)
endif()


## Add a configuration define for each library dependency (optional and built-in)  
## that we might need to know about in the code

set(DEPS_OPTIONS BOOST MPI OPENMP)             # vars of the form ENABLE_DEP
foreach(dep in ${DEPS_OPTIONS})
    if( ENABLE_${dep} )
        set(ATK_USE_${dep} TRUE)
    endif()
endforeach()

set(DEPS_BUILTIN CONDUIT HDF5 SPARSEHASH FMT)  # vars of the form DEP_FOUND
foreach(dep in ${DEPS_BUILTIN})
    if( ${dep}_FOUND )
        set(ATK_USE_${dep} TRUE)
    endif()
endforeach()


## Add a configuration define for each enabled toolkit component
set(COMPS COMMON LUMBERJACK SLIC SLAM SIDRE MINT QUEST SPIO)
foreach(comp in ${COMPS})
    if( ENABLE_${comp} )
        set(ATK_USE_${comp} TRUE)
    endif()
endforeach()




configure_file(
    components/common/src/config.hpp.in
    ${HEADER_INCLUDES_DIRECTORY}/common/config.hpp
)
