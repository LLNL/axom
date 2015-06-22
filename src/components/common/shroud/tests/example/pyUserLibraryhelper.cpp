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


PyObject *PP_ExClass1_to_Object(ExClass1 *grp)
{

    PyObject *voidobj;
    PyObject *args;
    PyObject *rv;

    voidobj = PyCapsule_New(grp, PY_ExClass1_capsule_name, NULL);
    args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, voidobj);
    rv = PyObject_Call((PyObject *) &PP_ExClass1_Type, args, NULL);
    Py_DECREF(args);
    return rv;
}

int PP_ExClass1_from_Object(PyObject *obj, void **addr)
{

    if (obj->ob_type != &PP_ExClass1_Type) {
	// raise exception
	return 0;	
    }
    PP_ExClass1 * self = (PP_ExClass1 *) obj;
    *addr = self->grp;
    
    return 1;

}

PyObject *PP_ExClass2_to_Object(ExClass2 *grp)
{

    PyObject *voidobj;
    PyObject *args;
    PyObject *rv;

    voidobj = PyCapsule_New(grp, PY_ExClass2_capsule_name, NULL);
    args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, voidobj);
    rv = PyObject_Call((PyObject *) &PP_ExClass2_Type, args, NULL);
    Py_DECREF(args);
    return rv;
}

int PP_ExClass2_from_Object(PyObject *obj, void **addr)
{

    if (obj->ob_type != &PP_ExClass2_Type) {
	// raise exception
	return 0;	
    }
    PP_ExClass2 * self = (PP_ExClass2 *) obj;
    *addr = self->grp;
    
    return 1;

}
}  // namespace example
}  // namespace nested
