! C code that will be inserted into Python module via shroud splicer blocks

// ----------------------------------------------------------------------
// ----- pySidremodule.hpp

// splicer begin include
#include "sidre/sidre.hpp"
// splicer end include

// ----------------------------------------------------------------------
// ----- pySidremodule.cpp


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


