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

###############################################################################
#
# file: src/examples/CMakeLists.txt
#
###############################################################################

# install using examples
install(DIRECTORY using-with-cmake using-with-make
        DESTINATION examples/axom)
# note: this is installed to examples/axom to avoid name conflicts with other
# libs (this will crop up if we use spack to install axom as a python 
# extension)
    

