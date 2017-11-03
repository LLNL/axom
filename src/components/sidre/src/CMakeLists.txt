#
# Specify all headers
#
set(sidre_headers
    core/Buffer.hpp
    core/Group.hpp
    core/DataStore.hpp
    core/View.hpp
    core/Attribute.hpp
    core/AttrValues.hpp
    core/Iterator.hpp
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
    core/AttrValues.cpp
    core/Iterator.cpp
    core/SidreUtilities.cpp
)

if(ENABLE_MPI)
   list(APPEND sidre_headers
        spio/IOManager.hpp
        spio/IOBaton.hpp
       )
   list(APPEND sidre_sources
        spio/IOManager.cpp
        spio/IOBaton.cpp
       )
endif()

if(SHROUD_FOUND)
    add_subdirectory(interface)
    if(ENABLE_MPI)
        add_subdirectory(spio/interface)
    endif()
endif()

if (ENABLE_FORTRAN)
    # The CAPI is only built to provide a Fortran interface.  If no Fortran,
    # disable CAPI.
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

    if(ENABLE_MPI)
        list(APPEND sidre_headers
             spio/interface/c_fortran/wrapIOManager.h
        )
        list(APPEND sidre_fortran_sources
             spio/interface/c_fortran/wrapfspio.f
             spio/interface/c_fortran/wrapIOManager.cpp
        )
    endif()
endif()

set(sidre_depends
    axom_utils
    conduit
    conduit_relay
    hdf5
    slic)

if(ENABLE_MPI)
   list(APPEND sidre_depends mpi)
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

if(SPARSEHASH_FOUND)
    blt_add_target_definitions(TO                 sidre
                               TARGET_DEFINITIONS ${SPARSEHASH_COMPILER_DEFS})
    target_include_directories(sidre PRIVATE
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/thirdparty/sparsehash>)
endif()

install(TARGETS              sidre
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include/sidre
        )
install(FILES       ${sidre_headers}
        DESTINATION include/sidre
        )
if(ENABLE_FORTRAN)
    install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/axom_sidre.mod DESTINATION lib/fortran)
endif()

# Only needed if building the Fortran interface
# XXX this can be used to automatically generate the files
#if (ENABLE_FORTRAN)
#    add_dependencies(copy_headers_sidre sidre_generate)
#endif()

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
    install(TARGETS              sidre
            EXPORT               ${PROJECT_NAME}-targets
            DESTINATION          lib
            INCLUDES DESTINATION include/sidre
            )
    install(FILES       ${sidre_headers}
            DESTINATION include/sidre
            )

#http://stackoverflow.com/questions/2152077/is-it-possible-to-get-cmake-to-build-both-a-static-and-shared-version-of-the-sam
# object library to avoid compile files twice
endif()
