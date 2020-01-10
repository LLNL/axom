# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

blt_add_executable(NAME lumberjack_basic_ex
                   SOURCES basicExample.cpp
                   OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                   DEPENDS_ON lumberjack 
                   FOLDER axom/lumberjack/examples )

