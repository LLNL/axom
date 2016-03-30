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

add_shroud( YAML_INPUT_FILE shroud_sidre.yaml
            C_FORTRAN_OUTPUT_DIR c_fortran
            PYTHON_OUTOUT_DIR python
            DEPENDS_SOURCE
                c_fortran/csidresplicer.c c_fortran/fsidresplicer.f
                python/pysidresplicer.c
            DEPENDS_BINARY genfsidresplicer.f
)

add_custom_command(
    OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.f
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/c_fortran/genfsidresplicer.py
    COMMAND ${PYTHON_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/c_fortran/genfsidresplicer.py fortran > 
            ${CMAKE_CURRENT_BINARY_DIR}/genfsidresplicer.f
)

# The CAPI is only built to provide a Fortran interface.  If no Fortran,
# disable CAPI.
if (ENABLE_FORTRAN)
    list(APPEND sidre_headers
        c_fortran/sidre.h
        c_fortran/SidreTypes.h
        # generated headers
        c_fortran/wrapSidre.h
        c_fortran/wrapDataStore.h
        c_fortran/wrapDataGroup.h
        c_fortran/wrapDataBuffer.h
        c_fortran/wrapDataView.h
    )   
    set(sidre_fortran_sources
        # generated source
        c_fortran/wrapSidre.cpp
        c_fortran/wrapDataStore.cpp
        c_fortran/wrapDataGroup.cpp
        c_fortran/wrapDataBuffer.cpp
        c_fortran/wrapDataView.cpp
        c_fortran/wrapfsidre.f
    )

    # wrappers depend on shroud runtime library
    set(c_fortran_depends_on
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
                  "${sidre_fortran_sources}"
              DEPENDS_ON
                  ${c_fortran_depends_on}
                  common
                  conduit
                  slic
             )
  
             
set_target_properties(sidre  PROPERTIES COMPILE_FLAGS "-Wno-unused-parameter")

# XXX - would like to have shroudrt loaded automatically with sidre
#target_link_libraries(sidre ${c_fortran_depends_on})

install(TARGETS sidre DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

# Only needed if building the Fortran interface
# XXX this can be used to automatically generate the files
#if (ENABLE_FORTRAN)
#    add_dependencies(copy_headers_sidre sidre_generate)
#endif()


if(ENABLE_PYTHON_API)
    set(sidre_python_sources
        python/pyDataStoretype.cpp
        python/pyDataGrouptype.cpp
        python/pyDataBuffertype.cpp
        python/pyDataViewtype.cpp
        python/pySidremodule.hpp
        python/pySidremodule.cpp
        python/pySidrehelper.cpp
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
