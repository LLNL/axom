// pystringsmodule.hpp
// This is generated code, do not edit
#ifndef PYSTRINGSMODULE_HPP
#define PYSTRINGSMODULE_HPP
#include "strings.hpp"
#include <Python.h>
#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif
// splicer begin header.include
// splicer end header.include
// splicer begin header.C_declaration
// splicer end header.C_declaration

// helper functions


extern PyObject *PY_error_obj;

#ifdef __cplusplus
extern "C" {
#endif
#ifdef IS_PY3K
#define MOD_INITBASIS PyInit_strings
#else
#define MOD_INITBASIS initstrings
#endif
PyMODINIT_FUNC MOD_INITBASIS(void);
#ifdef __cplusplus
}
#endif

#endif  /* PYSTRINGSMODULE_HPP */
