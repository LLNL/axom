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

blt_add_executable(NAME lumberjack_basic_ex
                   SOURCES basicExample.cpp
                   OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                   DEPENDS_ON lumberjack 
                   FOLDER axom/lumberjack/examples )

