// pyQUESTmodule.hpp
// This is generated code, do not edit
//
// Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
//
// Produced at the Lawrence Livermore National Laboratory
//
// LLNL-CODE-741217
//
// All rights reserved.
//
// This file is part of Axom.
//
// For details about use and distribution, please read axom/LICENSE.
//
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
