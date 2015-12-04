#
# makefile to build tutorial
#
WRAPDIR = ../tests/tutorial

FC = gfortran
FFLAGS = -g -ffree-form
CXXFLAGS = -g -I. -I$(WRAPDIR)

OBJS = \
    tutorial.o \
    wrapTutorial.o \
    wrapClass1.o \
    wrapftutorial.o \
    fstr_mod.o


VPATH = $(WRAPDIR):../runtime

all : testc testf

testc : testc.o libtutorial.a
	$(CXX) $^ -L. -ltutorial -o $@

testf : testf.o libtutorial.a
	$(FC) $^ -L. -ltutorial -lstdc++ -o $@

tutorial.o : tutorial.hpp

testf.o : wrapftutorial.o
wrapftutorial.o : fstr_mod.o

libtutorial.a : $(OBJS)
	ar cr $@ $?

.PHONY : clean
clean :
	rm -rf *~ *.o *.mod *.a testc testf
