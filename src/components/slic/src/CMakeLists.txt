#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

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
   if (ENABLE_lumberjack)
      list(APPEND slic_headers streams/LumberjackStream.hpp)
      list(APPEND slic_sources streams/LumberjackStream.cpp)
   endif()
endif()

# XXX ENABLE_FORTRAN_WRAPPER
if (ENABLE_FORTRAN)
    list(APPEND slic_headers
        interface/shroud_api/wrapSLIC.h
    )   
    set(slic_fortran_sources
        interface/shroud_api/wrapSLIC.cpp
        interface/shroud_api/shroudrt.cpp
        interface/shroud_api/wrapfslic.f
    )
    add_subdirectory(interface/shroud_api)

    # Must tell fortran where to look for modules
    include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
endif()

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( ${PROJECT_NAME}
                     "${slic_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )

#
# make the library
#
set(slic_depends_on
    common
)
if (ENABLE_lumberjack)
  list(APPEND slic_depends_on lumberjack)
endif()

make_library( NAME
                  slic

              SOURCES
                  "${slic_sources}"
                  "${slic_headers}"
                  "${slic_fortran_sources}"

              DEPENDS_ON
                  "${slic_depends_on}"
              )

install( TARGETS slic DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
