# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Examples for Lumberjack component
#------------------------------------------------------------------------------

blt_add_executable(NAME lumberjack_basic_ex
                   SOURCES basicExample.cpp
                   OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                   DEPENDS_ON axom 
                   FOLDER axom/lumberjack/examples )

