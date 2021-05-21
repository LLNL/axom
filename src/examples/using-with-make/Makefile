#------------------------------------------------------------------------------
# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
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
#------------------------------------------------------------------------------


INC_FLAGS=-I$(AXOM_DIR)/include/
LINK_FLAGS=-L$(AXOM_DIR)/lib/ -laxom

main:
	$(CXX) $(INC_FLAGS) example.cpp $(LINK_FLAGS) -o example


