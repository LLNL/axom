###############################################################################
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#
##############################################################################

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
    mpi
    )

if(SCR_FOUND)
     list(APPEND spio_depends scr)
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

if (ENABLE_FORTRAN)
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
install(TARGETS              spio
        EXPORT               ${PROJECT_NAME}-targets
        DESTINATION          lib
        INCLUDES DESTINATION include/spio
        )
install(FILES       ${spio_headers}
        DESTINATION include/spio
        )
if (ENABLE_FORTRAN)
    install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/axom_spio.mod DESTINATION lib/fortran)
endif()
