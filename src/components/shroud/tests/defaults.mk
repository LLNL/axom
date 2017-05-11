# shroud/tests/defaults.mk

compiler = gcc
#compiler = intel

# paths need the trailing slash
gcc.path = /usr/apps/gnu/4.7.1/bin/
gcc.path = /usr/apps/gnu/4.9.3/bin/
intel.path = /usr/local/tools/ic-15.0.187/bin/
intel.path = /usr/local/tools/ic-16.0.109/bin/

ifeq ($(compiler),gcc)
CC = $(gcc.path)gcc
CFLAGS = -g -Wall
CXX = $(gcc.path)g++
CXXFLAGS = -g -Wall
FC = $(gcc.path)gfortran
FFLAGS = -g -Wall -ffree-form
LIBS = -lstdc++
SHARED = -fPIC
LD_SHARED = -shared
endif

ifeq ($(compiler),intel)
CC = $(intel.path)icc
CFLAGS = -g
CXX = $(intel.path)icpc
CXXFLAGS = -g 
FC = $(intel.path)ifort
FFLAGS = -g -free
LIBS = -lstdc++
SHARED = -fPIC
LD_SHARED = -shared
endif

# 2.7
PYTHON_VER := $(shell $(PYTHON) -c "import sys;sys.stdout.write('{v[0]}.{v[1]}'.format(v=sys.version_info))")
PYTHON_PREFIX := $(shell $(PYTHON) -c "import sys;sys.stdout.write(sys.exec_prefix)")
PYTHON_BIN := $(PYTHON)
PYTHON_INC := -I$(PYTHON_PREFIX)/include/python$(PYTHON_VER)
PYTHON_LIB := -L$(PYTHON_PREFIX)/lib/python$(PYTHON_VER)/config -lpython$(PYTHON_VER) -ldl -lutil

LUA_PREFIX = $(abspath $(dir $(LUA))/..)
LUA_BIN = $(LUA)
LUA_INC = -I$(LUA_PREFIX)/include
LUA_LIB = -L$(LUA_PREFIX)/lib -llua -ldl

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $*.o $^

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $*.o $^

%.o %.mod  : %.f
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $*.o $^

%.o %.mod  : %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $*.o $^
