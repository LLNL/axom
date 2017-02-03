# shroud/tests/defaults.mk

CC = gcc
CFLAGS = -g -Wall
CXX = g++
CXXFLAGS = -g -Wall
FC = gfortran
FFLAGS = -g -Wall -ffree-form
LIBS = -lstdc++
SHARED = -fPIC
LD_SHARED = -shared

PYTHON = /home/taylor16/tpl/v2
PYTHON_BIN = $(PYTHON)/bin/python2
PYTHON_INC = -I$(PYTHON)/include/python2.7
PYTHON_LIBS = -L$(PYTHON)/lib/python2.7/config -lpython2.7 -ldl -lutil

LUA = /home/taylor16/tpl/v2
LUA_BIN = $(LUA)/bin/lua
LUA_INC = -I$(LUA)/include
LUA_LIBS = -L$(LUA)/lib -llua -ldl

%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $*.o $^

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $*.o $^

%.o %.mod  : %.f
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $*.o $^

%.o %.mod  : %.f90
	$(FC) $(FFLAGS) $(INCLUDE) -c -o $*.o $^
