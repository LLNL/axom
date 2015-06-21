! C code that will be inserted into Python module via shroud splicer blocks

// ----------------------------------------------------------------------
// ----- pySidremodule.hpp

// splicer begin include
#include "sidre/sidre.hpp"
// splicer end include

// splicer begin C_declaration
extern const char * datagroup_capsule_name;

extern PyObject *PP_DataGroup_to_Object(DataGroup *grp);

// splicer end C_declaration



// splicer begin class.DataStore.C_object
DataStore * ds;
// splicer end class.DataStore.C_object

// splicer begin class.DataGroup.C_object
DataGroup * grp;
// splicer end class.DataGroup.C_object

// splicer begin class.DataBuffer.C_object
DataBuffer * buf;
// splicer end class.DataBuffer.C_object

// splicer begin class.DataView.C_object
DataView * view;
// splicer end class.DataView.C_object



// ----------------------------------------------------------------------
// ----- pySidremodule.cpp

// splicer begin C_definition
const char * datagroup_capsule_name = "DataGroup";
// splicer end C_definition


// splicer begin additional_functions
PyObject *PP_DataGroup_to_Object(DataGroup *grp)
{
    PyObject *voidobj;
    PyObject *args;
    PyObject *rv;

    voidobj = PyCapsule_New(grp, datagroup_capsule_name, NULL);
    args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, voidobj);
    rv = PyObject_Call((PyObject *) &PY_DataGroup_Type, args, NULL);
    Py_DECREF(args);
    return rv;
}

// can be used with PyArg_ParseTupleAndKeywords().
// Pull DataGroup pointer out of Py_DataGroup object.
int PP_DataGroup_converter(PyObject *obj, void **addr)
{
    if (obj->ob_type != &PY_DataGroup_Type) {
	// raise exception
	return 0;	
    }
    PY_DataGroup * self = (PY_DataGroup *) obj;
    *addr = self->grp;
    
    return 1;
}

// splicer end additional_functions



// ----------------------------------------------------------------------
// ----- pyDataStoretype.cpp

// splicer begin class.DataStore.type.init
DataStore * ds = new DataStore();
self->ds = ds;
return 0;
// splicer end class.DataStore.type.init

// splicer begin class.DataStore.method.getRoot
DataGroup * grp = self->ds->getRoot();
PyObject *rv = PP_DataGroup_to_Object(grp);
return rv;
// splicer end class.DataStore.method.getRoot

// ----------------------------------------------------------------------
// ----- pyDataGrouptype.cpp

// splicer begin class.DataGroup.type.init
    PyObject *grpobj;

    /* By requiring a PyCapsule, it is difficult to call directly from Python.
     * But the C++ constructors are private so that makes sense.
     */
    if (!PyArg_ParseTuple(args, "O!:DataGroup_init",
                          &PyCapsule_Type, &grpobj))
        return -1;

    /* capsule_dbnode */
    DataGroup *grp = static_cast<DataGroup *>(PyCapsule_GetPointer(grpobj, datagroup_capsule_name));
    self->grp = grp;
    if (grp == NULL && PyErr_Occurred())
	return -1;
    
    return 0;
// splicer end class.DataGroup.type.init


// splicer begin class.DataGroup.method.getName
const std::string & name = self->grp->getName();
PyObject * rv = PyString_FromString(name.c_str());
return rv;
// splicer end class.DataGroup.method.getName


// ----------------------------------------------------------------------
// ----- pyDataBuffertype.cpp


// ----------------------------------------------------------------------
// ----- pyDataViewtype.cpp


