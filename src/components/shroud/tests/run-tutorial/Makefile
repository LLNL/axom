#
# run-tutorial
#
include $(top)/tests/defaults.mk

INCLUDE = \
    -I$(top)/tests/tutorial \
    -I$(top)/tests/run-tutorial \
    -I$(top)/tests/fruit

VPATH = \
    $(top)/tests/tutorial \
    $(top)/tests/run-tutorial \
    $(top)/tests/fruit

OBJS = \
    tutorial.o \
    wrapftutorial.o \
    wrapTutorial.o \
    wrapClass1.o \
    shroudrt.o \
    fruit.o \
    main.o

tutorial : $(OBJS)
	$(FC) $(FFLAGS) $(LIBS) $^ -o $@
