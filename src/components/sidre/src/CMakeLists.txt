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

# Must tell fortran where to look for modules
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

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
if (ENABLE_FORTRAN)
    add_dependencies(copy_headers_sidre sidre_generate)
endif()
