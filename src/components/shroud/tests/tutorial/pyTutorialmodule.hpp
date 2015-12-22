// pyTutorialmodule.hpp
// This is generated code, do not edit
#ifndef PYTUTORIALMODULE_HPP
#define PYTUTORIALMODULE_HPP
#include "tutorial.hpp"
#include <Python.h>
#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif
// splicer begin header.include
// splicer end header.include
namespace tutorial {
extern PyTypeObject PY_Class1_Type;
// splicer begin header.C_declaration
// splicer end header.C_declaration

// helper functions
extern const char *PY_Class1_capsule_name;
PyObject *PP_Class1_to_Object(Class1 *addr);
int PP_Class1_from_Object(PyObject *obj, void **addr);

// splicer begin class.Class1.C_declaration
// splicer end class.Class1.C_declaration

typedef struct {
PyObject_HEAD
    Class1 * BBB;
    // splicer begin class.Class1.C_object
    // splicer end class.Class1.C_object
} PY_Class1;

extern PyObject *PY_error_obj;

#ifdef __cplusplus
extern "C" {
#endif
#ifdef IS_PY3K
#define MOD_INITBASIS PyInit_tutorial
#else
#define MOD_INITBASIS inittutorial
#endif
PyMODINIT_FUNC MOD_INITBASIS(void);
#ifdef __cplusplus
}
#endif

}  // namespace tutorial
#endif  /* PYTUTORIALMODULE_HPP */
