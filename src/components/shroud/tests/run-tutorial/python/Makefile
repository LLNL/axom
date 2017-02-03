#
#  Compile the python tutorial
#
include $(top)/tests/defaults.mk

INCLUDE = \
    -I$(top)/tests/tutorial \
    -I$(top)/tests/run-tutorial \
    $(PYTHON_INC)

VPATH = \
    $(top)/tests/tutorial \
    $(top)/tests/run-tutorial \
    $(top)/tests/run-tutorial/python

OBJS = \
    tutorial.o \
    pyTutorialmodule.o \
    pyClass1type.o \
    pyTutorialhelper.o

CXXFLAGS += $(SHARED)

all : tutorial.so simple

tutorial.so : $(OBJS)
	$(CXX) $(LD_SHARED) -o $@ $^ $(LIBS)

simple : testpython.o $(OBJS)
	$(CXX) -pthread -o $@ $^ $(PYTHON_LIBS)
#	g++ -pthread -o $@ $^ $(PYTHON_LIBS)

clean :
	rm -f *.so *.o simple
.PHONY : clean

