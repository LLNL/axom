#
# Specify all headers
#
set(sidre_headers
    core/Buffer.hpp
    core/Group.hpp
    core/DataStore.hpp
    core/View.hpp
    core/Attribute.hpp
    core/MapCollection.hpp
    core/SidreTypes.hpp
    core/SidreDataTypeIds.h
    core/sidre.hpp
)

#
# Specify all sources
#
set(sidre_sources
    core/Buffer.cpp
    core/Group.cpp
    core/DataStore.cpp
    core/View.cpp
    core/Attribute.cpp
    core/SidreUtilities.cpp
)

if(SHROUD_FOUND)
    add_subdirectory(interface)
endif()

# The CAPI is only built to provide a Fortran interface.  If no Fortran,
# disable CAPI.
if (ENABLE_FORTRAN)
    list(APPEND sidre_headers
        interface/sidre.h
        interface/SidreTypes.h
        # generated headers
        interface/c_fortran/wrapSidre.h
        interface/c_fortran/wrapDataStore.h
        interface/c_fortran/wrapGroup.h
        interface/c_fortran/wrapBuffer.h
        interface/c_fortran/wrapView.h
    )   
    set(sidre_fortran_sources
        # generated source
        interface/c_fortran/wrapSidre.cpp
        interface/c_fortran/wrapDataStore.cpp
        interface/c_fortran/wrapGroup.cpp
        interface/c_fortran/wrapBuffer.cpp
        interface/c_fortran/wrapView.cpp
        interface/c_fortran/wrapfsidre.F
        interface/c_fortran/shroudrt.cpp
    )
endif()

set(sidre_depends
    axom_utils
    conduit
    conduit_relay
    hdf5
    slic)

if(SPARSEHASH_FOUND)
   set(sidre_depends ${sidre_depends} sparsehash)
endif()

#
# make the library
#
blt_add_library( NAME
                     sidre
                 SOURCES
                     "${sidre_sources}"
                     "${sidre_fortran_sources}"
                 HEADERS
                     "${sidre_headers}"
                 HEADERS_OUTPUT_SUBDIR
                     sidre
                 DEPENDS_ON
                     ${sidre_depends}
                 )

# Only needed if building the Fortran interface
# XXX this can be used to automatically generate the files
#if (ENABLE_FORTRAN)
#    add_dependencies(copy_headers_sidre sidre_generate)
#endif()

if(SPARSEHASH_FOUND)
   blt_add_target_definitions(TO
                              sidre
                              TARGET_DEFINITIONS
                              ${SPARSEHASH_COMPILER_DEFS})
endif()

if(FALSE)  # The Python module is not create yet.
#if(ENABLE_PYTHON_API)
    set(sidre_python_sources
        interface/python/pyDataStoretype.cpp
        interface/python/pyGrouptype.cpp
        interface/python/pyBuffertype.cpp
        interface/python/pyViewtype.cpp
        interface/python/pySidremodule.cpp
        interface/python/pySidrehelper.cpp
    )
    set(sidre_python_headers
        interface/python/pySidremodule.hpp
    )

    blt_add_library(
        NAME sidre
        PYTHON_MODULE
        SOURCES
            ${sidre_sources}
            ${sidre_python_sources}
        HEADERS
            ${sidre_headers}
            ${sidre_python_headers}
        DEPENDS_ON
            common
            conduit
            conduit_io
            hdf5
            slic
    )

#http://stackoverflow.com/questions/2152077/is-it-possible-to-get-cmake-to-build-both-a-static-and-shared-version-of-the-sam
# object library to avoid compile files twice
endif()
