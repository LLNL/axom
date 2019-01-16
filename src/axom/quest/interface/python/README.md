

To build the Python interface for quest

```
% cd $CMAKE_BINARY_DIR/axom/quest/interface/python
% python setup.py build
% setenv PYTHONPATH `pwd`/build/lib.linux-x86_64-2.7
% python
>>> import quest

% python3 setup.py build
% setenv PYTHONPATH `pwd`/build/lib.linux-x86_64-3.6
% python3
>>> import quest
```

`extra_compile_args` in `setpy.py` must include flags to enable C++11.

Some functions in `quest_shroud.yaml` are not wrapped pending fixes to shroud.
(version 0.10.1 is currently used)

Need to Axom build with shared libaries.

```
# Make shared lib
set(BUILD_SHARED_LIBS "ON" CACHE BOOL "")
```

Need location of MPI

```
set(MPI_HOME "/usr/tce/packages/mvapich2/mvapich2-2.2-intel-18.0.2" CACHE PATH "")
```


Building as static libraries is not working.  Need to list every library
but still results in a runtime error.  Loading with a static MPI library
would probably result in problems with mpi4py module.

```
    extra_link_args=["-L" + os.path.join(builddir, "lib"),
                     "-lquest", "-lmint", "-lslic", "-llumberjack"]
ImportError: /g/g14/taylor/axom/work/build-quest-only-debug/axom/quest/interface/python/build/lib.linux-x86_64-2.7/quest.so: undefined symbol: MPI_Barrier
                    
```


```
# Sidre requires conduit and hdf5, so disable it in this host-config
set(AXOM_ENABLE_SIDRE OFF CACHE BOOL "")


ImportError: /g/g14/taylor/axom/work/build-rzgenie-toss_3_x86_64_ib-clang@4.0.0-debug/axom/quest/interface/python/build/lib.linux-x86_64-2.7/quest.so: undefined symbol: _ZN4axom5sidre5Group16s_path_delimiterE
```
