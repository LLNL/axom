// pyQUESTmodule.hpp
// This is generated code, do not edit
//
// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef PYQUESTMODULE_HPP
#define PYQUESTMODULE_HPP
#include <Python.h>
// splicer begin header.include
// splicer end header.include

// splicer begin header.C_declaration
// splicer end header.C_declaration

// helper functions


extern PyObject* PY_error_obj;

#if PY_MAJOR_VERSION >= 3
extern "C" PyMODINIT_FUNC PyInit_quest(void);
#else
extern "C" PyMODINIT_FUNC initquest(void);
#endif

#endif  /* PYQUESTMODULE_HPP */
