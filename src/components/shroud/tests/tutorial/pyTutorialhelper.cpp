// pyTutorialhelper.cpp
// This is generated code, do not edit
#include "pyTutorialmodule.hpp"
namespace tutorial {
const char *PY_Class1_capsule_name = "Class1";


PyObject *PP_Class1_to_Object(Class1 *addr)
{
    // splicer begin class.Class1.helper.to_object
    PyObject *voidobj;
    PyObject *args;
    PyObject *rv;

    voidobj = PyCapsule_New(addr, PY_Class1_capsule_name, NULL);
    args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, voidobj);
    rv = PyObject_Call((PyObject *) &PY_Class1_Type, args, NULL);
    Py_DECREF(args);
    return rv;
    // splicer end class.Class1.helper.to_object
}

int PP_Class1_from_Object(PyObject *obj, void **addr)
{
    // splicer begin class.Class1.helper.from_object
    if (obj->ob_type != &PY_Class1_Type) {
        // raise exception
        return 0;	
    }
    PY_Class1 * self = (PY_Class1 *) obj;
    *addr = self->BBB;
    return 1;
    // splicer end class.Class1.helper.from_object
}
}  // namespace tutorial
