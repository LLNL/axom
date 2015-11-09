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
// splicer begin header.include
// splicer end header.include
namespace example {
namespace nested {
extern PyTypeObject PP_ExClass1_Type;
extern PyTypeObject PP_ExClass2_Type;
// splicer begin header.C_declaration
// splicer end header.C_declaration

// helper functions
extern const char *PY_ExClass1_capsule_name;
extern const char *PY_ExClass2_capsule_name;
PyObject *PP_ExClass1_to_Object(ExClass1 *addr);
int PP_ExClass1_from_Object(PyObject *obj, void **addr);
PyObject *PP_ExClass2_to_Object(ExClass2 *addr);
int PP_ExClass2_from_Object(PyObject *obj, void **addr);

// splicer begin class.ExClass1.C_declaration
// splicer end class.ExClass1.C_declaration

typedef struct {
PyObject_HEAD
    ExClass1 * BBB;
    // splicer begin class.ExClass1.C_object
    // splicer end class.ExClass1.C_object
} PP_ExClass1;
// splicer begin class.ExClass2.C_declaration
// splicer end class.ExClass2.C_declaration

typedef struct {
PyObject_HEAD
    ExClass2 * BBB;
    // splicer begin class.ExClass2.C_object
    // splicer end class.ExClass2.C_object
} PP_ExClass2;

extern PyObject *PP_error_obj;

#ifdef __cplusplus
extern "C" {
#endif
#ifdef IS_PY3K
#define MOD_INITBASIS PyInit_userlibrary
#else
#define MOD_INITBASIS inituserlibrary
#endif
PyMODINIT_FUNC MOD_INITBASIS(void);
#endif
#ifdef __cplusplus
}
#endif

}  // namespace example
}  // namespace nested
