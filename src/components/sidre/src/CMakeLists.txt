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
)

# The CAPI is only built to provide a Fortran interface.  If no Fortran,
# disable CAPI.
if (ENABLE_FORTRAN)
    list(APPEND sidre_headers
        capi/sidre.h
        capi/SidreTypes.h
        capi/SidreWrapperHelpers.hpp
        capi/SidreAllocatable.hpp
        # generated headers
        capi/wrapSidre.h
        capi/wrapDataStore.h
        capi/wrapDataGroup.h
        capi/wrapDataBuffer.h
        capi/wrapDataView.h
    )   
    list(APPEND sidre_sources
        # generated source
        capi/wrapSidre.cpp
        capi/wrapDataStore.cpp
        capi/wrapDataGroup.cpp
        capi/wrapDataBuffer.cpp
        capi/wrapDataView.cpp
        capi/wrapfsidre.f
        capi/SidreAllocatable.cpp
        capi/SidreAllocatablef.f
    )
    add_subdirectory(capi)

    # wrappers depend on shroud runtime library
    set(capi_depends_on
        shroudrt
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
                  "${sidre_headers}"
              DEPENDS_ON
                  ${capi_depends_on}
                  common
                  conduit
                  slic
             )
set_target_properties(sidre  PROPERTIES COMPILE_FLAGS "-Wno-unused-parameter")

# XXX - would like to have shroudrt loaded automatically with sidre
#target_link_libraries(sidre ${capi_depends_on})

install(TARGETS sidre DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

# Only needed if building the Fortran interface
# XXX this can be used to automatically generate the files
#if (ENABLE_FORTRAN)
#    add_dependencies(copy_headers_sidre sidre_generate)
#endif()
