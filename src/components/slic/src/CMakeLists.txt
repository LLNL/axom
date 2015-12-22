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
    make_library( LIBRARY_NAME
                      slic

                  LIBRARY_SOURCES
                      "${slic_sources}"
                      "${slic_headers}"
                  DEPENDS_ON
                      common
                  )


install( TARGETS slic DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
