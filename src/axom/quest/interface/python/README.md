

# Python interface for Quest signed distance

Build Axom as usual, then use the generated `setup.py` script.
It will work with Python 2.7 or 3.4+.

```
% cd $CMAKE_BINARY_DIR
% make
% cd axom/quest/interface/python

% python setup.py build
% setenv PYTHONPATH `pwd`/build/lib.linux-x86_64-2.7
% python quest_type.py

% python3 setup.py build
% setenv PYTHONPATH `pwd`/build/lib.linux-x86_64-3.6
% python3 quest_type.py
```

The file `quest_type.py` requires the module `mpi4py`.
Both `mpi4py` and Axom should be build with the same version of MPI.

`distutils` uses the same compiler and flags that were used to build Python.
`extra_compile_args` in `setpy.py` must include flags to enable C++11.
These flags may be different than the ones used to compile Axom.

Some functions in `quest_shroud.yaml` are not wrapped pending fixes to shroud.
(version 0.10.1 is currently used)

A sample host file is provided in `host-configs/quest-only.cmake`.
Need to build Axom with shared libaries and the location of MPI.

```
# Make shared lib
set(BUILD_SHARED_LIBS "ON" CACHE BOOL "")

set(MPI_HOME "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-18.0.2" CACHE PATH "")
```
