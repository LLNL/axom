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
(version 0.13.0 is currently used)

A sample host-config file is provided in
`host-configs/other/no-tpl-toss_3_x86_64_ib-intel@19.0.4.cmake`. This host-config requires
no third-party libraries but does not enable Sidre.  The Python interface requires you build
Axom with shared libraries which can be enabled via the CMake command line with
`-DBUILD_SHARED_LIBS=ON` and defining MPI_HOME, which is already enabled in the sample host-config.

To build only quest, and its dependencies:

    ./config-build.py -hc host-configs/other/no-tpl-toss_3_x86_64_ib-intel@19.0.4.cmake -DBUILD_SHARED_LIBS=ON -DAXOM_ENABLE_ALL_COMPONENTS=BOOL:OFF -DAXOM_ENABLE_QUEST=BOOL:ON  -DAXOM_ENABLE_MINT=BOOL:ON  -DAXOM_ENABLE_SLIC=BOOL:ON -DAXOM_ENABLE_SLAM=BOOL:ON  -DAXOM_ENABLE_PRIMAL=BOOL:ON  -DAXOM_ENABLE_SPIN=BOOL:ON

