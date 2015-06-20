#
# Compile Python module
#

PYTHON_INC = -I/home/taylor16/testzzz-ets-4.6.3-4.3f/include/python2.7
PYTHON_LIB = -L/home/taylor16/testzzz-ets-4.6.3-4.3f/lib/python2.7/config -lpython2.7.a

CXX = /home/taylor16/gapps/gcc-4.9.0/bin/g++
CXXFLAGS = -g -fPIC $(PYTHON_INC)





OBJS = \
    python_module.o \
    python_datastore.o \
    python_datagroup.o \
    python_databuffer.o \
    python_dataview.o


sidre.so : $(OBJS)
	gcc -shared -fPIC -Wl,-soname,$@ -o $@ $(OBJS)
#	$(LD) $(OBJS) --export-dynamic -o $@


clean : 
	rm -f *.o *.so
.PHONY : clean


