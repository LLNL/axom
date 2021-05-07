[comment]: # (#################################################################)
[comment]: # (Copyright 2017-2021, Lawrence Livermore National Security, LLC)
[comment]: # (and Axom Project Developers. See the top-level LICENSE file)
[comment]: # (for details.)
[comment]: #
[comment]: # (# SPDX-License-Identifier: BSD-3-Clause)
[comment]: # (#################################################################)


This directory contains some batch scripts to build third-party libraries 
or source for multiple configurations on different clusters at LLNL.

About
-----
There are two types of batch scripts for each network and cluster type:
* The `src` scripts call `<axom>/scripts/llnl_scripts/build_src.py` 
  to build and test all configurations on the given cluster.
* The `all_compilers` scripts call `<axom>/scripts/llnl_scripts/build_tpls.py` 
  to build Axom's third-party libraries (via uberenv) on the given cluster.
	As part of this process, it build and test the code using each of these configurations .

Usage
-----
To use these scripts, log into a cluster on an LLNL network (`CZ` or `RZ`) and call the appropriate script:
* On `toss3` (e.g. `ruby`) use the `msub` scripts: `msub msub_llnl_{c,r}z_toss3_{src,all_compilers}.sh`
* On `blueos` (e.g. `lassen`) use the `bsub` scripts: `bsub < bsub_llnl_{c,r}z_blueos_{src,all_compilers}.sh`
