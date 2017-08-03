#
# Specify all headers
#
set(spio_headers
    IOBaton.hpp
    IOManager.hpp
    )

set(spio_depends
    hdf5
    conduit
    sidre
    fmt
    )

if(SCR_FOUND)
     set(spio_depends
         ${spio_depends}
         scr
         )
endif()

#
# Specify all sources
#
set(spio_sources
    IOBaton.cpp
    IOManager.cpp
    )

if(SHROUD_FOUND)
    add_subdirectory(interface)
endif()

if(ENABLE_FORTRAN)
    set(spio_fortran_headers
        # generated headers
        interface/c_fortran/wrapIOManager.h
    )
    set(spio_fortran_sources
        # generated source
        interface/c_fortran/wrapIOManager.cpp
        interface/c_fortran/wrapfspio.f
        interface/c_fortran/shroudrt.cpp
    )
endif()

#
# make the library
#
blt_add_library( NAME
                     spio

                 SOURCES
                     "${spio_sources}"
                     "${spio_fortran_sources}"

                 HEADERS
                     "${spio_headers}"
                     "${spio_fortran_headers}"

                 HEADERS_OUTPUT_SUBDIR
                     spio

                 DEPENDS_ON
                     ${spio_depends}
                )
