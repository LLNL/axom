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
    Logger.hpp
    LogStream.hpp
    MessageLevel.h
    UnitTestLogger.hpp
    slic.hpp

    logstreams/GenericOutputStream.hpp
    )

#
# Specify all sources
#
set(slic_sources
    Logger.cpp
    LogStream.cpp

    logstreams/GenericOutputStream.cpp
    )

if (ENABLE_MPI)
   list(APPEND slic_headers logstreams/SynchronizedStream.hpp)
   list(APPEND slic_sources logstreams/SynchronizedStream.cpp)
endif()

#
# make the library
#
if ( ENABLE_MPI )

    make_library( LIBRARY_NAME
                      slic
                  LIBRARY_SOURCES
                      "${slic_sources}"
                      "${slic_headers}"
                  WITH_MPI
                  )

else()

    make_library( LIBRARY_NAME
                      slic
                  LIBRARY_SOURCES
                      "${slic_sources}"
                      "${slic_headers}"
                  )

endif()
add_dependencies(slic copy_headers_common)
add_dependencies(slic copy_headers_slic)

install( TARGETS slic DESTINATION lib EXPORT ${PROJECT_NAME}-targets)

#
# custom target to copy header files in the "include" folder of the build tree.
#
copy_headers_target( ${PROJECT_NAME}
                     "${slic_headers}"
                     ${HEADER_INCLUDES_DIRECTORY}/${PROJECT_NAME} )
