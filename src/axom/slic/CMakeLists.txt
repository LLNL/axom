#------------------------------------------------------------------------------
# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

################################
# Specify all headers/sources
################################
set(slic_headers
    core/Logger.hpp
    core/LogStream.hpp
    core/MessageLevel.h
    core/UnitTestLogger.hpp

    interface/slic.hpp

    internal/stacktrace.hpp

    streams/GenericOutputStream.hpp
    )

set(slic_sources
    core/Logger.cpp
    core/LogStream.cpp

    interface/slic.cpp

    internal/stacktrace.cpp

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
        interface/c_fortran/typesSLIC.h
    )
    list(APPEND slic_sources
        interface/c_fortran/wrapSLIC.cpp
        interface/c_fortran/wrapfslic.f
    )
endif()

################################
# Make/Install the library
################################
set(slic_depends_on core)

if (AXOM_ENABLE_LUMBERJACK)
    list(APPEND slic_depends_on lumberjack mpi)
endif()

if (WIN32)
    list(APPEND slic_depends_on dbghelp)  
endif()

blt_add_library(NAME        slic
                SOURCES     ${slic_sources}
                HEADERS     ${slic_headers}
                DEPENDS_ON  ${slic_depends_on}
                FOLDER      axom/slic
                OBJECT      TRUE
                )

axom_write_unified_header(NAME    slic
                          HEADERS ${slic_headers})

axom_install_component(NAME    slic
                       HEADERS ${slic_headers}
                       )


################################
# Add examples
################################
if (AXOM_ENABLE_EXAMPLES)
    add_subdirectory(examples)
endif()

################################
# Add tests
################################
if (AXOM_ENABLE_TESTS)
    add_subdirectory(tests)
endif()

################################
# Add docs
################################
if (AXOM_ENABLE_DOCS)
    add_subdirectory(docs)
endif()

################################
# Add code checks
################################
axom_add_code_checks(PREFIX slic)
