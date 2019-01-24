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

################################
# documentation targets
################################

blt_add_sphinx_target( mint_docs )

if ( DOXYGEN_FOUND )
  add_dependencies( mint_docs doxygen_docs )
endif()
