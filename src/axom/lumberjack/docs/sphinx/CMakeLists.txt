# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

blt_add_sphinx_target(lumberjack_docs)

if(DOXYGEN_FOUND)
   add_dependencies(lumberjack_docs doxygen_docs)
endif()
