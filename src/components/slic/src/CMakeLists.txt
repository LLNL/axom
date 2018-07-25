#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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
#------------------------------------------------------------------------------

#
# Specify all headers
#
set(slic_headers
    core/Logger.hpp
    core/LogStream.hpp
    core/MessageLevel.h
    core/UnitTestLogger.hpp

    interface/slic.hpp

    streams/GenericOutputStream.hpp
    )

#
# Specify all sources
#
set(slic_sources
    core/Logger.cpp
    core/LogStream.cpp

    interface/slic.cpp

    streams/GenericOutputStream.cpp
    )

if (ENABLE_MPI)
    list(APPEND slic_headers streams/SynchronizedStream.hpp)
    list(APPEND slic_sources streams/SynchronizedStream.cpp)
    if (AXOM_ENABLE_LUMBERJACK)
        list(APPEND slic_headers streams/LumberjackStream.hpp)
        list(APPEND slic_sources streams/LumberjackStream.cpp)
    endif()
endif()

if(SHROUD_FOUND)
    add_subdirectory(interface)
endif()

if (ENABLE_FORTRAN)
    list(APPEND slic_headers
        interface/c_fortran/wrapSLIC.h
    )
    set(slic_fortran_sources
        interface/c_fortran/wrapSLIC.cpp
        interface/c_fortran/wrapfslic.f
    )
endif()

#
# make the library
#
set(slic_depends_on
    axom_utils
)
if (AXOM_ENABLE_LUMBERJACK)
    list(APPEND slic_depends_on lumberjack mpi)
endif()

if (WIN32)
    list(APPEND slic_depends_on dbghelp)  
endif()

blt_add_library( NAME
                    slic

                 SOURCES
                    "${slic_sources}"
                    "${slic_fortran_sources}"

                 HEADERS
                    "${slic_headers}"

                 HEADERS_OUTPUT_SUBDIR
                     slic

                 DEPENDS_ON
                    "${slic_depends_on}"

                 FOLDER
                     axom/slic
              )
install(TARGETS              slic
        EXPORT               slic-targets
        DESTINATION          lib
        INCLUDES DESTINATION include/slic
        )
install(FILES       ${slic_headers}
        DESTINATION include/slic
        )
if(ENABLE_FORTRAN)
    install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/axom_slic.mod DESTINATION lib/fortran)
endif()
