#
#  Compile the lua tutorial
#
include $(top)/tests/defaults.mk

INCLUDE = \
    -I$(top)/tests/tutorial \
    -I$(top)/tests/run-tutorial \
    $(LUA_INC)

VPATH = \
    $(top)/tests/tutorial \
    $(top)/tests/run-tutorial \
    $(top)/tests/run-tutorial/lua

OBJS = \
    tutorial.o \
    luaTutorialmodule.o

CXXFLAGS += $(SHARED)

all : tutorial.so # simple

tutorial.so : $(OBJS)
	$(CXX) $(LD_SHARED) -o $@ $^ $(LIBS)

# XXX testlua.c has some 5.2 features
simple : testlua.o $(OBJS)
	$(CC) -o $@ $^ $(LUA_LIB)

clean :
	rm -f *.so *.o simple
.PHONY : clean

print-debug:
	@echo LUA=$(LUA)
	@echo LUA_PREFIX=$(LUA_PREFIX)
	@echo LUA_INC=$(LUA_INC)
	@echo LUA_LIB=$(LUA_LIB)
