
WRAPDIR = ../tests/tutorial

FC = gfortran
FFLAGS = -g -ffree-form
CXXFLAGS = -g -I. -I$(WRAPDIR)

OBJS = \
    tutorial.o \
    wrapTutorial.o \
    wrapftutorial.o \
    fstr_mod.o


VPATH = $(WRAPDIR)

ftest : ftest.o libtutorial.a
	$(FC) $^ -L. -ltutorial -o $@

ftest.o : wrapftutorial.o
wrapftutorial.o : fstr_mod.o

libtutorial.a : $(OBJS)
	ar cr $@ $?

.PHONY : clean
clean :
	rm -rf *~ *.o *.mod *.a ftest
