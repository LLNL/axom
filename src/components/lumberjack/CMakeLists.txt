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

################################
# Specify all headers/sources
################################
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

set(lumberjack_sources
    BinaryTreeCommunicator.cpp
    Lumberjack.cpp
    Message.cpp
    MPIUtility.cpp
    RootCommunicator.cpp
    )


################################
# Make/Install the library
################################
blt_add_library( NAME
                     lumberjack
                 SOURCES
                     ${lumberjack_sources}
                 HEADERS
                     ${lumberjack_headers}
                 DEPENDS_ON
                     axom_utils mpi
                 FOLDER
                     axom/lumberjack
                 )

install(TARGETS              lumberjack
        EXPORT               lumberjack-targets
        DESTINATION          lib
        INCLUDES DESTINATION include/lumberjack
        )
install(FILES       ${lumberjack_headers}
        DESTINATION include/lumberjack
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
axom_add_code_checks(PREFIX lumberjack)

################################
# Create CMake importable
# exports for all of our targets
################################
install(EXPORT lumberjack-targets DESTINATION lib/cmake)
