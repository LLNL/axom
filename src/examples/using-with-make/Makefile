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
# Example that shows how to use an installed instance of Axom in Makefile
# based build system.
#
# To build:
#  env AXOM_DIR={axom install path} make
#  ./example
#
# Note:
#  If you use features that leverage external third party libraries 
#  (ex: hdf5 features), you may need to pass additional include 
#   and linking flags. 
# 
#
###############################################################################


INC_FLAGS=-I$(AXOM_DIR)/include/
LINK_FLAGS=-L$(AXOM_DIR)/lib/ -laxom

main:
	$(CXX) $(INC_FLAGS) example.cpp $(LINK_FLAGS) -o example


