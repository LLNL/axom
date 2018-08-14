## How to update third party libraries

```
cd scripts/llnl_scripts/batch_scripts
make clean    # remove old log scripts
make tpl-rz 
```

A log file for each batch job will be created.
After the jobs are finished search for `SUCCESS` or `ERROR`.

The third-party-libraries will be build in the directory
`/usr/workspace/wsrzc/axom/thirdparty_libs/builds/YYYY_MM_DD_HH_HH_MM`.
When `uberenv` finished it will create some files with the suffix
`.cmake` which must be copied into the `host-configs` directory.

Then repeat for the CZ.  The libraries will be built in
`/usr/workspace/wsa/axom/thirdparty_libs/builds/YYYY_MM_DD_HH_HH_MM`.

## How to test source for all compilers

```
cd scripts/llnl_scripts/batch_scripts
make clean    # remove old log scripts
make src-rz 
```

The source will be built and tested in the root of the cloned
directory with the name `_axom_build_and_test_YYYY_MM_DD_HH_HH_MM`.
