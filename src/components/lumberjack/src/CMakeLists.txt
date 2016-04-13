#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

configure_file(config.hpp.in "${HEADER_INCLUDES_DIRECTORY}/lumberjack/config.hpp")

#
# Specify all headers
#
set(lumberjack_headers
    BinaryTreeCommunicator.hpp
    Combiner.hpp
    Communicator.hpp
    Lumberjack.hpp
    Message.hpp
    MPIUtility.hpp
    RootCommunicator.hpp
    TextEqualityCombiner.hpp
    )

#
# Specify all sources
#
set(lumberjack_sources
    BinaryTreeCommunicator.cpp
    Lumberjack.cpp
    Message.cpp
    MPIUtility.cpp
    RootCommunicator.cpp
    )


#
# make the library
#
blt_add_library( NAME
                     lumberjack
                 SOURCES
                     "${lumberjack_sources}"
                 HEADERS
                     "${lumberjack_headers}"
                 DEPENDS_ON
                     common
                 )

install(TARGETS lumberjack DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
