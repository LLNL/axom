// pySidremodule.hpp
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
// splicer begin include
#include "sidre/sidre.hpp"
// splicer end include
namespace asctoolkit {
namespace sidre {
extern PyTypeObject PY_DataStore_Type;
extern PyTypeObject PY_DataGroup_Type;
extern PyTypeObject PY_DataBuffer_Type;
extern PyTypeObject PY_DataView_Type;
// splicer begin C_declaration
extern const char * datagroup_capsule_name;

extern PyObject *PP_DataGroup_to_Object(DataGroup *grp);

// splicer end C_declaration
// splicer begin class.DataStore.C_declaration
// splicer end class.DataStore.C_declaration

typedef struct {
PyObject_HEAD
// splicer begin class.DataStore.C_object
DataStore * ds;
// splicer end class.DataStore.C_object
} PY_DataStore;
// splicer begin class.DataGroup.C_declaration
// splicer end class.DataGroup.C_declaration

typedef struct {
PyObject_HEAD
// splicer begin class.DataGroup.C_object
DataGroup * grp;
// splicer end class.DataGroup.C_object
} PY_DataGroup;
// splicer begin class.DataBuffer.C_declaration
// splicer end class.DataBuffer.C_declaration

typedef struct {
PyObject_HEAD
// splicer begin class.DataBuffer.C_object
DataBuffer * buf;
// splicer end class.DataBuffer.C_object
} PY_DataBuffer;
// splicer begin class.DataView.C_declaration
// splicer end class.DataView.C_declaration

typedef struct {
PyObject_HEAD
// splicer begin class.DataView.C_object
DataView * view;
// splicer end class.DataView.C_object
} PY_DataView;

extern PyObject *PY_error_obj;

#ifdef __cplusplus
extern "C" {
#endif
#ifdef IS_PY3K
#define MOD_INITBASIS PyInit_sidre
#else
#define MOD_INITBASIS initsidre
#endif
PyMODINIT_FUNC MOD_INITBASIS(void);
#endif
#ifdef __cplusplus
}
#endif

}  // namespace asctoolkit
}  // namespace sidre
