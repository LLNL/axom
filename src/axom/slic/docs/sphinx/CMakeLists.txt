#------------------------------------------------------------------------------
# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------

################################
# documentation targets
################################

# add a target to generate documentation with Sphinx

blt_add_sphinx_target( slic_docs )

if (DOXYGEN_FOUND)
  add_dependencies(slic_docs doxygen_docs)
endif()
