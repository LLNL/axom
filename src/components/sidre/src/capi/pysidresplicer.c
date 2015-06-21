! C code that will be inserted into Python module via shroud splicer blocks

// splicer begin include
#include "sidre/sidre.hpp"
// splicer end include


Global variables
// splicer begin C_definition
const char * datagroup_capsule_name = "DataGroup";
// splicer end C_definition


// splicer begin C_declaration
extern const char * datagroup_capsule_name;
// splicer end C_declaration


// splicer begin class.DataStore.C_object
DataStore * ds;
// splicer end class.DataStore.C_object

// splicer begin class.DataStore.type.init
DataStore * ds = new DataStore();
self->ds = ds;
return 0;
// splicer end class.DataStore.type.init

// splicer begin class.DataStore.method.getRoot
DataGroup * grp = self->ds->getRoot();
PyObject *voidobj;
PyObject *args0;
PyObject *rv;

voidobj = PyCapsule_New(grp, datagroup_capsule_name, NULL);
args0 = PyTuple_New(1);
PyTuple_SET_ITEM(args, 0, voidobj);
rv = PyObject_Call((PyObject *) &PY_DataGroup_Type, args0, NULL);
Py_DECREF(args0);
return rv;
// splicer end class.DataStore.method.getRoot


    BA_dbnode *node;
    PyObject *nodeobj;

    /* By requiring a PyCapsule, it is difficult to call directly
     * from Python. Intended to be called from PB_dbnode_to_instance.
     */
    if (!PyArg_ParseTuple(args, "O!:dbnode_init",
                          &PyCapsule_Type, &nodeobj))
        return -1;

    /* capsule_dbnode */
    node = (BA_dbnode *) PyCapsule_GetPointer(nodeobj, "BA_dbnode"); 
    self->node = node;
    if (node == NULL && PyErr_Occurred())
	return -1;
    
    return 0;
