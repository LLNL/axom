###############################################################################
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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
###############################################################################

#
# Compile Python module
#
CMAKE_BINARY_DIR = /home/taylor16/datastore/sprint6/build-debug

PYTHON_INC = -I/home/taylor16/testzzz-ets-4.6.3-4.3f/include/python2.7
PYTHON_LIB = -L/home/taylor16/testzzz-ets-4.6.3-4.3f/lib/python2.7/config -lpython2.7.a

CXX = /home/taylor16/gapps/gcc-4.9.0/bin/g++
CXXFLAGS = -g -fPIC $(PYTHON_INC) -I$(CMAKE_BINARY_DIR)/include

LDFLAGS = -L$(CMAKE_BINARY_DIR)/lib -lsidre -lslic -lconduit -lstdc++
LDFLAGS = -L$(CMAKE_BINARY_DIR)/lib -Wl,--whole-archive -lsidre -lslic -lconduit -Wl,--no-whole-arc -lstdc++ -lgfortran



OBJS = \
	pySidremodule.o \
	pyDataStoretype.o \
	pyDataGrouptype.o \
	pyDataBuffertype.o \
	pyDataViewtype.o \
	pySidrehelper.o

all : clean generate sidre.so


sidre.so : $(OBJS)
	gcc -shared -fPIC -Wl,-soname,$@ -o $@ $(LDFLAGS) $(OBJS)
#	$(LD) $(OBJS) --export-dynamic -o $@


generate :
	( touch api.yaml; cd $(CMAKE_BINARY_DIR); make sidre_generate )

clean : 
	rm -f *.o *.so
.PHONY : clean


