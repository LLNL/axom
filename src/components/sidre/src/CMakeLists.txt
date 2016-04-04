#
# Specify all headers
#
add_subdirectory(interface)

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

    # Must tell fortran where to look for modules
    include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
endif()
#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( ${PROJECT_NAME} "${sidre_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )

#
# make the library
#
make_library( LIBRARY_NAME
                  sidre
              LIBRARY_SOURCES
                  "${sidre_sources}"
                  "${sidre_fortran_sources}"
              DEPENDS_ON
                  common
                  conduit
                  conduit_io
                  slic
             )
  
             
set_target_properties(sidre  PROPERTIES COMPILE_FLAGS "-Wno-unused-parameter")

install(TARGETS sidre DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

# Only needed if building the Fortran interface
# XXX this can be used to automatically generate the files
#if (ENABLE_FORTRAN)
#    add_dependencies(copy_headers_sidre sidre_generate)
#endif()


if(ENABLE_PYTHON_API)
    set(sidre_python_sources
        interface/python/pyDataStoretype.cpp
        interface/python/pyDataGrouptype.cpp
        interface/python/pyDataBuffertype.cpp
        interface/python/pyDataViewtype.cpp
        interface/python/pySidremodule.hpp
        interface/python/pySidremodule.cpp
        interface/python/pySidrehelper.cpp
    )

    include_directories(${PYTHON_INCLUDE_DIRS})

    set(BUILD_SHARED_LIBS TRUE)
    make_library(
        LIBRARY_NAME
            sidre2
        LIBRARY_SOURCES
            ${sidre_sources}
            ${sidre_python_sources}
        DEPENDS_ON
            common
            conduit
            slic
    )

    # Python wants the name without leading 'lib'
    set_target_properties(sidre2
        PROPERTIES
          PREFIX ""
    )
    set_target_properties(sidre2 PROPERTIES OUTPUT_NAME sidre)

#http://stackoverflow.com/questions/2152077/is-it-possible-to-get-cmake-to-build-both-a-static-and-shared-version-of-the-sam
# object library to avoid compile files twice
endif(ENABLE_PYTHON_API)
