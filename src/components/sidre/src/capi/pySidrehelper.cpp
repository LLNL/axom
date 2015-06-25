// pySidrehelper.cpp
// This is generated code, do not edit
//
// Copyright (c) 2015, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory.
//
// All rights reserved.
//
// This source code cannot be distributed without permission and
// further review from Lawrence Livermore National Laboratory.
//
#include "pySidremodule.hpp"
namespace asctoolkit {
namespace sidre {
const char *PY_DataStore_capsule_name = "DataStore";
const char *PY_DataGroup_capsule_name = "DataGroup";
const char *PY_DataBuffer_capsule_name = "DataBuffer";
const char *PY_DataView_capsule_name = "DataView";


PyObject *PP_DataStore_to_Object(DataStore *addr)
{
    // splicer begin class.DataStore.helper.to_object
    PyObject *voidobj;
    PyObject *args;
    PyObject *rv;
    
    voidobj = PyCapsule_New(addr, PY_DataStore_capsule_name, NULL);
    args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, voidobj);
    rv = PyObject_Call((PyObject *) &PY_DataStore_Type, args, NULL);
    Py_DECREF(args);
    return rv;
    // splicer end class.DataStore.helper.to_object
}

int PP_DataStore_from_Object(PyObject *obj, void **addr)
{
    // splicer begin class.DataStore.helper.from_object
    if (obj->ob_type != &PY_DataStore_Type) {
        // raise exception
        return 0;	
    }
    PY_DataStore * self = (PY_DataStore *) obj;
    *addr = self->BBB;
    return 1;
    // splicer end class.DataStore.helper.from_object
}

PyObject *PP_DataGroup_to_Object(DataGroup *addr)
{
    // splicer begin class.DataGroup.helper.to_object
    PyObject *voidobj;
    PyObject *args;
    PyObject *rv;
    
    voidobj = PyCapsule_New(addr, PY_DataGroup_capsule_name, NULL);
    args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, voidobj);
    rv = PyObject_Call((PyObject *) &PY_DataGroup_Type, args, NULL);
    Py_DECREF(args);
    return rv;
    // splicer end class.DataGroup.helper.to_object
}

int PP_DataGroup_from_Object(PyObject *obj, void **addr)
{
    // splicer begin class.DataGroup.helper.from_object
    if (obj->ob_type != &PY_DataGroup_Type) {
        // raise exception
        return 0;	
    }
    PY_DataGroup * self = (PY_DataGroup *) obj;
    *addr = self->BBB;
    return 1;
    // splicer end class.DataGroup.helper.from_object
}

PyObject *PP_DataBuffer_to_Object(DataBuffer *addr)
{
    // splicer begin class.DataBuffer.helper.to_object
    PyObject *voidobj;
    PyObject *args;
    PyObject *rv;
    
    voidobj = PyCapsule_New(addr, PY_DataBuffer_capsule_name, NULL);
    args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, voidobj);
    rv = PyObject_Call((PyObject *) &PY_DataBuffer_Type, args, NULL);
    Py_DECREF(args);
    return rv;
    // splicer end class.DataBuffer.helper.to_object
}

int PP_DataBuffer_from_Object(PyObject *obj, void **addr)
{
    // splicer begin class.DataBuffer.helper.from_object
    if (obj->ob_type != &PY_DataBuffer_Type) {
        // raise exception
        return 0;	
    }
    PY_DataBuffer * self = (PY_DataBuffer *) obj;
    *addr = self->BBB;
    return 1;
    // splicer end class.DataBuffer.helper.from_object
}

PyObject *PP_DataView_to_Object(DataView *addr)
{
    // splicer begin class.DataView.helper.to_object
    PyObject *voidobj;
    PyObject *args;
    PyObject *rv;
    
    voidobj = PyCapsule_New(addr, PY_DataView_capsule_name, NULL);
    args = PyTuple_New(1);
    PyTuple_SET_ITEM(args, 0, voidobj);
    rv = PyObject_Call((PyObject *) &PY_DataView_Type, args, NULL);
    Py_DECREF(args);
    return rv;
    // splicer end class.DataView.helper.to_object
}

int PP_DataView_from_Object(PyObject *obj, void **addr)
{
    // splicer begin class.DataView.helper.from_object
    if (obj->ob_type != &PY_DataView_Type) {
        // raise exception
        return 0;	
    }
    PY_DataView * self = (PY_DataView *) obj;
    *addr = self->BBB;
    return 1;
    // splicer end class.DataView.helper.from_object
}
}  // namespace asctoolkit
}  // namespace sidre
