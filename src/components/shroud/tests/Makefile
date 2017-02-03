#
# testsdir passed in from parent makefile
#
FC = gfortran
FC_FREE = -ffree-form
LIBS = -lstdc++

#
# MAKECMDGOALS will be tutorial or strings
#
INCLUDE = \
    -I$(testsdir)/$(MAKECMDGOALS) \
    -I$(testsdir)/run-$(MAKECMDGOALS) \
    -I$(testsdir)/fruit

VPATH = \
    $(testsdir)/$(MAKECMDGOALS) \
    $(testsdir)/run-$(MAKECMDGOALS) \
    $(testsdir)/fruit

tutorial_objs = \
    tutorial.o \
    wrapftutorial.o \
    wrapTutorial.o \
    wrapClass1.o \
    shroudrt.o \
    fruit.o \
    main.o

strings_objs = \
    strings.o \
    wrapfstrings.o \
    wrapstrings.o \
    fruit.o \
    main.o

all : tutorial strings


#tutorial: VPATH = $(testsdir)/tutorial $(testsdir)/run-tutorial $(testsdir)/fruit
tutorial : $(tutorial_objs)
	$(FC) $(FFLAGS) $(LIBS) $^ -o $@

#strings : VPATH = $(testsdir)/strings $(testsdir)/run-strings $(testsdir)/fruit
strings : $(strings_objs)
	echo $(VPATH)
	$(FC) $(FFLAGS) $(LIBS) $^ -o $@

clean :
	rm -rf *.o *.mod tutorial strings

.PHONY : all clean


%.o : %.c
	$(CC) $(CPPFLAGS) $(INCLUDE) -c -o $*.o $^

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $*.o $^

%.o %.mod  : %.f
	$(FC) $(FC_FREE) $(FFLAGS) -c -o $*.o $^

%.o %.mod  : %.f90
	$(FC) $(FC_FREE) $(FFLAGS) -c -o $*.o $^
