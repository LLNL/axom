// pySLICmodule.hpp
// This is generated code, do not edit
//
// Copyright (c) 2015, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
//
// All rights reserved.
//
// This source code cannot be distributed without permission and
// further review from Lawrence Livermore National Laboratory.
//

/*
 * This is generated code.
 * Any edits must be made between the splicer.begin and splicer.end blocks.
 * All other edits will be lost.
 * Once a block is edited remove the 'UNMODIFIED' on the splicer.begin
 * comment to allow the block to be preserved when it is regenerated.
 */

#ifndef HDR_BASISMODULE
#define HDR_BASISMODULE
#include <Python.h>
#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif
// splicer begin header.include
// splicer end header.include
namespace asctoolkit {
namespace slic {
// splicer begin header.C_declaration
// splicer end header.C_declaration

// helper functions


extern PyObject *PY_error_obj;

#ifdef __cplusplus
extern "C" {
#endif
#ifdef IS_PY3K
#define MOD_INITBASIS PyInit_slic
#else
#define MOD_INITBASIS initslic
#endif
PyMODINIT_FUNC MOD_INITBASIS(void);
#endif
#ifdef __cplusplus
}
#endif

}  // namespace asctoolkit
}  // namespace slic
