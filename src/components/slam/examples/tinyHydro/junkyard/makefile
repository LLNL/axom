PYTHON_DIR=/usr/gapps/asctoolkit/thirdparty_libs/spack/opt/chaos_5_x86_64_ib/gcc@4.4.7/python@2.7.8
PYTHON_EXE=$(PYTHON_DIR)/bin/python2.7

CC=g++
INCLUDE = -I$(PYTHON_DIR)/include/python2.7
# Optimized with profile
# CFLAGS = -fPIC -O3 $(INCLUDE)
# # Optimized
CFLAGS = -fPIC -O3  $(INCLUDE)
# Debug
# CFLAGS = -g -fPIC  $(INCLUDE)
#PYBINDGEN = /g/g16/doug/python/PyBindGen-0.17.0 # not actually using this
LIBS =  -lrt -l$(PYTHON_DIR)/lib/libpython2.7.so

obs: vectorXY.so PolygonMeshXY.so State.so Hydro.so Part.so

mini: main.o PolygonMeshXY.o State.o HydroC.o 
	$(CC) -o $@ $^ -lrt

wc: VectorXY.hpp PolygonMeshXY.[ch] State.[ch] HydroC.[ch] vectorXY-binding.py PolygonMeshXY-binding.py State-binding.py HydroC-binding.py Part.[ch] Part-binding.py
	wc $^

###### level 3 ################## (requires level 2 )

Hydro.so: HydroC-binding.o HydroC.o State.o PolygonMeshXY.o Part.o
	$(CC) -shared -o $@ $^ $(LIBS)

HydroC-binding.cpp: HydroC-binding.py HydroC.hpp


###### level 2 ################## (requires level 1)

State.so: State-binding.o State.o  VectorXY.hpp PolygonMeshXY.hpp Part.o 
	$(CC) -shared -o $@ $^ $(LIBS)

State-binding.cpp: State-binding.py State.hpp

###### level 1 ################## (requires level 0)

PolygonMeshXY.so: PolygonMeshXY-binding.o PolygonMeshXY.o VectorXY.hpp
	$(CC) -shared -o $@ $^ $(LIBS)

PolygonMeshXY-binding.cpp: PolygonMeshXY-binding.py PolygonMeshXY.hpp

###### level 0 ##################  (requires nothing)

Part.so: Part-binding.o Part.o 
	$(CC) -shared -o $@ $^ $(LIBS)

Part-binding.cpp: Part-binding.py Part.hpp

vectorXY.so: vectorXY-binding.o
	$(CC) -shared $^ -o $@

vectorXY-binding.cpp: vectorXY-binding.py VectorXY.hpp
	$(PYTHON_EXE) $< > $@

%.o : %.cpp
	$(CC) -c $(CFLAGS) $< -o $@


%.cpp: %.py
	$(PYTHON_EXE) $< > $@

clean:
	rm *.o *.so *-binding.cpp

# reminders:
# $^ = "all the dependents"
# $< = "first dependent"
# $@ = "target"
