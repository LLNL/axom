// pytestnameshelper.cpp
// This is generated code, do not edit
#include "pytestnamesmodule.hpp"
const char *PY_Names_capsule_name = "Names";


PyObject *PP_Names_to_Object(Names *addr)
{
    // splicer begin class.Names.helper.to_object
    PyObject *voidobj;
    PyObject *args;
    PyObject *rv;

    voidobj = PyCapsule_New(addr, PY_Names_capsule_name, NULL);
    args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, voidobj);
    rv = PyObject_Call((PyObject *) &PY_Names_Type, args, NULL);
    Py_DECREF(args);
    return rv;
    // splicer end class.Names.helper.to_object
}

int PP_Names_from_Object(PyObject *obj, void **addr)
{
    // splicer begin class.Names.helper.from_object
    if (obj->ob_type != &PY_Names_Type) {
        // raise exception
        return 0;
    }
    PY_Names * self = (PY_Names *) obj;
    *addr = self->BBB;
    return 1;
    // splicer end class.Names.helper.from_object
}
