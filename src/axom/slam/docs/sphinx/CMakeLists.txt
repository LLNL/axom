# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Sphinx targets for slam
#------------------------------------------------------------------------------

blt_add_sphinx_target( slam_docs )

if (DOXYGEN_FOUND)
   add_dependencies(slam_docs doxygen_docs)
endif()
