#
# Specify all headers
#

# The CAPI is only built to provide a Fortran interface.  If no Fortran,
# disable CAPI.
if (ENABLE_FORTRAN)
    add_subdirectory(capi)
endif()

set(sidre_headers
    core/Collections.hpp
    core/DataBuffer.hpp
    core/DataGroup.hpp
    core/DataStore.hpp
    core/DataView.hpp
    core/SidreTypes.hpp
    core/sidre.hpp
    ${capi_headers}
)

#
# Specify all sources
#
set(sidre_sources
    core/DataBuffer.cpp
    core/DataGroup.cpp
    core/DataStore.cpp
    core/DataView.cpp
    ${capi_source}
)

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
                common
                conduit
                slic
             )

set_target_properties(sidre  PROPERTIES Fortran_FORMAT "FREE")

install(TARGETS sidre DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

# Only needed if building the Fortran interface
if (ENABLE_FORTRAN)
    add_dependencies(copy_headers_sidre sidre_generate)
endif()
