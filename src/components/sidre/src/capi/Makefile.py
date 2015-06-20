#
# Compile Python module
#

PYTHON_INC = -I/home/taylor16/testzzz-ets-4.6.3-4.3f/include/python2.7
PYTHON_LIB = -L/home/taylor16/testzzz-ets-4.6.3-4.3f/lib/python2.7/config -lpython2.7.a

CXX = /home/taylor16/gapps/gcc-4.9.0/bin/g++
CXXFLAGS = -g -fPIC $(PYTHON_INC) -I/home/taylor16/datastore/sprint4/build-debug/include


OBJS = \
	pySidremodule.o \
	pyDataStoretype.o \
	pyDataGrouptype.o \
	pyDataBuffertype.o \
	pyDataViewtype.o

sidre.so : $(OBJS)
	gcc -shared -fPIC -Wl,-soname,$@ -o $@ $(OBJS)
#	$(LD) $(OBJS) --export-dynamic -o $@


clean : 
	rm -f *.o *.so
.PHONY : clean


