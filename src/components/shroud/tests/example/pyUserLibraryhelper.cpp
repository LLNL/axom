// pyUserLibraryhelper.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
#include "pyUserLibrarymodule.hpp"
namespace example {
namespace nested {
const char *PY_ExClass1_capsule_name = "ExClass1";
const char *PY_ExClass2_capsule_name = "ExClass2";


PyObject *PP_ExClass1_to_Object(ExClass1 *addr)
{
    // splicer begin class.ExClass1.helper.to_object
    PyObject *voidobj;
    PyObject *args;
    PyObject *rv;

    voidobj = PyCapsule_New(addr, PY_ExClass1_capsule_name, NULL);
    args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, voidobj);
    rv = PyObject_Call((PyObject *) &PP_ExClass1_Type, args, NULL);
    Py_DECREF(args);
    return rv;
    // splicer end class.ExClass1.helper.to_object
}

int PP_ExClass1_from_Object(PyObject *obj, void **addr)
{
    // splicer begin class.ExClass1.helper.from_object
    if (obj->ob_type != &PP_ExClass1_Type) {
        // raise exception
        return 0;
    }
    PP_ExClass1 * self = (PP_ExClass1 *) obj;
    *addr = self->BBB;
    return 1;
    // splicer end class.ExClass1.helper.from_object
}

PyObject *PP_ExClass2_to_Object(ExClass2 *addr)
{
    // splicer begin class.ExClass2.helper.to_object
    PyObject *voidobj;
    PyObject *args;
    PyObject *rv;

    voidobj = PyCapsule_New(addr, PY_ExClass2_capsule_name, NULL);
    args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, voidobj);
    rv = PyObject_Call((PyObject *) &PP_ExClass2_Type, args, NULL);
    Py_DECREF(args);
    return rv;
    // splicer end class.ExClass2.helper.to_object
}

int PP_ExClass2_from_Object(PyObject *obj, void **addr)
{
    // splicer begin class.ExClass2.helper.from_object
    if (obj->ob_type != &PP_ExClass2_Type) {
        // raise exception
        return 0;
    }
    PP_ExClass2 * self = (PP_ExClass2 *) obj;
    *addr = self->BBB;
    return 1;
    // splicer end class.ExClass2.helper.from_object
}
}  // namespace example
}  // namespace nested
