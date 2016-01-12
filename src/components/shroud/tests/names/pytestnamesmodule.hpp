// pytestnamesmodule.hpp
// This is generated code, do not edit
#ifndef PYTESTNAMESMODULE_HPP
#define PYTESTNAMESMODULE_HPP
#include <Python.h>
#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif
// splicer begin header.include
// splicer end header.include
extern PyTypeObject PY_Names_Type;
// splicer begin header.C_declaration
// splicer end header.C_declaration

// helper functions
extern const char *PY_Names_capsule_name;
PyObject *PP_Names_to_Object(Names *addr);
int PP_Names_from_Object(PyObject *obj, void **addr);

// splicer begin class.Names.C_declaration
// splicer end class.Names.C_declaration

typedef struct {
PyObject_HEAD
    Names * BBB;
    // splicer begin class.Names.C_object
    // splicer end class.Names.C_object
} PY_Names;

extern PyObject *PY_error_obj;

#ifdef __cplusplus
extern "C" {
#endif
#ifdef IS_PY3K
#define MOD_INITBASIS PyInit_testnames
#else
#define MOD_INITBASIS inittestnames
#endif
PyMODINIT_FUNC MOD_INITBASIS(void);
#ifdef __cplusplus
}
#endif

#endif  /* PYTESTNAMESMODULE_HPP */
