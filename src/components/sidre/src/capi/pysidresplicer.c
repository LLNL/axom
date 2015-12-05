! C code that will be inserted into Python module via shroud splicer blocks

// ----------------------------------------------------------------------
// ----- pySidremodule.hpp

// splicer begin header.include
#include "sidre/sidre.hpp"
#include "sidre/SidreTypes.h"
#include "SidreWrapperHelpers.hpp"
// splicer end header.include

// ----------------------------------------------------------------------
// ----- pySidremodule.cpp

// splicer begin C_init_body
PyModule_AddIntConstant(m, "InvalidIndex", -1);
// splicer end C_init_body

// ----------------------------------------------------------------------
// ----- pyDataStoretype.cpp

// splicer begin class.DataStore.type.init
DataStore * ds = new DataStore();
self->BBB = ds;
return 0;
// splicer end class.DataStore.type.init

// splicer begin class.DataStore.type.richcompare
PyObject *rv = Py_NotImplemented;
if (PyObject_IsInstance(other, (PyObject*) &PY_DataStore_Type)) {
    PY_DataStore *pyother = (PY_DataStore *) other;
    switch (opid) {
    case Py_EQ:
	if (self->BBB == pyother->BBB) {
	    rv = Py_True;
	} else {
	    rv = Py_False;
	}
	break;
    case Py_NE:
	if (self->BBB != pyother->BBB) {
	    rv = Py_True;
	} else {
	    rv = Py_False;
	}
	break;
    case Py_LT:
    case Py_LE:
    case Py_GE:
    case Py_GT:
	break;
    }
 }
Py_INCREF(rv);
return rv;
// splicer end class.DataStore.type.richcompare
}

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
    DataGroup *grp = static_cast<DataGroup *>(PyCapsule_GetPointer(grpobj, PY_DataGroup_capsule_name));
    self->BBB = grp;
    if (grp == NULL && PyErr_Occurred())
	return -1;
    
    return 0;
// splicer end class.DataGroup.type.init


// splicerX begin class.DataGroup.method.getName
const std::string & name = self->BBB->getName();
PyObject * rv = PyString_FromString(name.c_str());
return rv;
// splicerX end class.DataGroup.method.getName

// splicer begin class.DataGroup.type.richcompare
PyObject *rv = Py_NotImplemented;
if (PyObject_IsInstance(other, (PyObject*) &PY_DataGroup_Type)) {
    PY_DataGroup *pyother = (PY_DataGroup *) other;
    switch (opid) {
    case Py_EQ:
	if (self->BBB == pyother->BBB) {
	    rv = Py_True;
	} else {
	    rv = Py_False;
	}
	break;
    case Py_NE:
	if (self->BBB != pyother->BBB) {
	    rv = Py_True;
	} else {
	    rv = Py_False;
	}
	break;
    case Py_LT:
    case Py_LE:
    case Py_GE:
    case Py_GT:
	break;
    }
 }
Py_INCREF(rv);
return rv;
// splicer end class.DataGroup.type.richcompare

// ----------------------------------------------------------------------
// ----- pyDataBuffertype.cpp

// splicer begin class.DataBuffer.type.richcompare
PyObject *rv = Py_NotImplemented;
if (PyObject_IsInstance(other, (PyObject*) &PY_DataBuffer_Type)) {
    PY_DataBuffer *pyother = (PY_DataBuffer *) other;
    switch (opid) {
    case Py_EQ:
	if (self->BBB == pyother->BBB) {
	    rv = Py_True;
	} else {
	    rv = Py_False;
	}
	break;
    case Py_NE:
	if (self->BBB != pyother->BBB) {
	    rv = Py_True;
	} else {
	    rv = Py_False;
	}
	break;
    case Py_LT:
    case Py_LE:
    case Py_GE:
    case Py_GT:
	break;
    }
 }
Py_INCREF(rv);
return rv;
// splicer end class.DataBuffer.type.richcompare

// ----------------------------------------------------------------------
// ----- pyDataViewtype.cpp

// splicer begin class.DataView.type.richcompare
PyObject *rv = Py_NotImplemented;
if (PyObject_IsInstance(other, (PyObject*) &PY_DataView_Type)) {
    PY_DataView *pyother = (PY_DataView *) other;
    switch (opid) {
    case Py_EQ:
	if (self->BBB == pyother->BBB) {
	    rv = Py_True;
	} else {
	    rv = Py_False;
	}
	break;
    case Py_NE:
	if (self->BBB != pyother->BBB) {
	    rv = Py_True;
	} else {
	    rv = Py_False;
	}
	break;
    case Py_LT:
    case Py_LE:
    case Py_GE:
    case Py_GT:
	break;
    }
 }
Py_INCREF(rv);
return rv;
// splicer end class.DataView.type.richcompare
