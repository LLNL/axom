#
# run-strings
#
include $(top)/tests/defaults.mk

INCLUDE = \
    -I$(top)/tests/strings \
    -I$(top)/tests/run-strings \
    -I$(top)/tests/fruit

VPATH = \
    $(top)/tests/strings \
    $(top)/tests/run-strings \
    $(top)/tests/fruit

OBJS = \
    strings.o \
    wrapfstrings.o \
    wrapstrings.o \
    shroudrt.o \
    fruit.o \
    main.o

strings : $(OBJS)
	$(FC) $(FFLAGS) $(LIBS) $^ -o $@
