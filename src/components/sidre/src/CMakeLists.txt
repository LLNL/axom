#
# Specify all headers
#
set(sidre_headers
    core/Collections.hpp
    core/DataBuffer.hpp
    core/DataGroup.hpp
    core/DataStore.hpp
    core/DataView.hpp
    core/SidreTypes.hpp
    core/SidreDataTypeIds.h
    core/sidre.hpp
)

#
# Specify all sources
#
set(sidre_sources
    core/DataBuffer.cpp
    core/DataGroup.cpp
    core/DataStore.cpp
    core/DataView.cpp
    core/SidreUtilities.cpp
)

if(ENABLE_SHROUD)
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
        interface/c_fortran/wrapDataGroup.h
        interface/c_fortran/wrapDataBuffer.h
        interface/c_fortran/wrapDataView.h
    )   
    set(sidre_fortran_sources
        # generated source
        interface/c_fortran/wrapSidre.cpp
        interface/c_fortran/wrapDataStore.cpp
        interface/c_fortran/wrapDataGroup.cpp
        interface/c_fortran/wrapDataBuffer.cpp
        interface/c_fortran/wrapDataView.cpp
        interface/c_fortran/wrapfsidre.f
        interface/c_fortran/shroudrt.cpp
    )
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
                     common
                     conduit
                     conduit_relay
                     hdf5
                     slic
                 )
  
             
set_target_properties(sidre  PROPERTIES COMPILE_FLAGS ${ATK_DISABLE_UNUSED_PARAMETER_WARNINGS})

# Only needed if building the Fortran interface
# XXX this can be used to automatically generate the files
#if (ENABLE_FORTRAN)
#    add_dependencies(copy_headers_sidre sidre_generate)
#endif()


if(FALSE)  # The Python module is not create yet.
#if(ENABLE_PYTHON_API)
    set(sidre_python_sources
        interface/python/pyDataStoretype.cpp
        interface/python/pyDataGrouptype.cpp
        interface/python/pyDataBuffertype.cpp
        interface/python/pyDataViewtype.cpp
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
