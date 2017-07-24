#-------------------------------------------------------------------------------
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

################################
# documentation targets
################################

# add a target to generate documentation with Sphinx

add_sphinx_target( primal_docs )

if (DOXYGEN_FOUND)
  add_dependencies(primal_docs primal_doxygen)
endif()
