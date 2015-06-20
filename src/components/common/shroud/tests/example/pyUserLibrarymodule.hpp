// pyUserLibrarymodule.hpp
// This is generated code, do not edit
// blah blah
// yada yada
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
// splicer begin C_declaration
// splicer end C_declaration
// splicer begin class.C_declaration
// splicer end class.C_declaration

typedef struct {
PyObject_HEAD
// splicer begin class.C_object
// splicer end class.C_object
} PP_ExClass1;
// splicer begin class.C_declaration
// splicer end class.C_declaration

typedef struct {
PyObject_HEAD
// splicer begin class.C_object
// splicer end class.C_object
} PP_ExClass2;

extern PyObject *PP_error_obj;

#ifdef IS_PY3K
#define MOD_INITBASIS PyInit_userlibrary
#else
#define MOD_INITBASIS inituserlibrary
#endif
PyMODINIT_FUNC MOD_INITBASIS(void);
#endif

