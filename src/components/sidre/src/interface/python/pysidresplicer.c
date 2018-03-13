!C code that will be inserted into Python module via shroud splicer blocks

/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


// ----------------------------------------------------------------------
// ----- pySidremodule.hpp

// splicer begin header.include
#include "sidre/sidre.hpp"
#include "sidre/SidreTypes.h"
// splicer end header.include

// ----------------------------------------------------------------------
// ----- pySidremodule.cpp

// splicer begin C_init_body
PyModule_AddIntConstant(m, "InvalidIndex", -1);

PyModule_AddIntConstant(m, "NO_TYPE_ID", SIDRE_NO_TYPE_ID);
PyModule_AddIntConstant(m, "INT8_ID ", SIDRE_INT8_ID);
PyModule_AddIntConstant(m, "INT16_ID", SIDRE_INT16_ID);
PyModule_AddIntConstant(m, "INT32_ID", SIDRE_INT32_ID);
PyModule_AddIntConstant(m, "INT64_ID", SIDRE_INT64_ID);

PyModule_AddIntConstant(m, "UINT8_ID ", SIDRE_UINT8_ID);
PyModule_AddIntConstant(m, "UINT16_ID", SIDRE_UINT16_ID);
PyModule_AddIntConstant(m, "UINT32_ID", SIDRE_UINT32_ID);
PyModule_AddIntConstant(m, "UINT64_ID", SIDRE_UINT64_ID);

PyModule_AddIntConstant(m, "FLOAT32_ID  ", SIDRE_FLOAT32_ID);
PyModule_AddIntConstant(m, "FLOAT64_ID  ", SIDRE_FLOAT64_ID);
PyModule_AddIntConstant(m, "CHAR8_STR_ID", SIDRE_CHAR8_STR_ID);

PyModule_AddIntConstant(m, "INT_ID", SIDRE_INT_ID);
PyModule_AddIntConstant(m, "UINT_ID", SIDRE_UINT_ID);
PyModule_AddIntConstant(m, "LONG_ID", SIDRE_LONG_ID);
PyModule_AddIntConstant(m, "ULONG_ID", SIDRE_ULONG_ID);
PyModule_AddIntConstant(m, "FLOAT_ID", SIDRE_FLOAT_ID);
PyModule_AddIntConstant(m, "DOUBLE_ID", SIDRE_DOUBLE_ID);
// splicer end C_init_body

// ----------------------------------------------------------------------
// ----- pyDataStoretype.cpp

// splicer begin class.DataStore.type.init
DataStore* ds = new DataStore();
self->BBB = ds;
return 0;
// splicer end class.DataStore.type.init

// splicer begin class.DataStore.type.richcompare
PyObject* rv = Py_NotImplemented;
if (PyObject_IsInstance(other, (PyObject*) &PY_DataStore_Type))
{
  PY_DataStore* pyother = (PY_DataStore*) other;
  switch (opid)
  {
  case Py_EQ:
    if (self->BBB == pyother->BBB)
    {
      rv = Py_True;
    }
    else
    {
      rv = Py_False;
    }
    break;
  case Py_NE:
    if (self->BBB != pyother->BBB)
    {
      rv = Py_True;
    }
    else
    {
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
PyObject* grpobj;

/* By requiring a PyCapsule, it is difficult to call directly from Python.
 * But the C++ constructors are private so that makes sense.
 */
if (!PyArg_ParseTuple(args, "O!:DataGroup_init",
                      &PyCapsule_Type, &grpobj))
  return -1;

/* capsule_dbnode */
DataGroup* grp = static_cast<DataGroup *>
                 (PyCapsule_GetPointer(grpobj, PY_DataGroup_capsule_name));
self->BBB = grp;
if (grp == NULL && PyErr_Occurred())
  return -1;

return 0;
// splicer end class.DataGroup.type.init


// splicerX begin class.DataGroup.method.getName
const std::string & name = self->BBB->getName();
PyObject* rv = PyString_FromString(name.c_str());
return rv;
// splicerX end class.DataGroup.method.getName

// splicer begin class.DataGroup.type.richcompare
PyObject* rv = Py_NotImplemented;
if (PyObject_IsInstance(other, (PyObject*) &PY_DataGroup_Type))
{
  PY_DataGroup* pyother = (PY_DataGroup*) other;
  switch (opid)
  {
  case Py_EQ:
    if (self->BBB == pyother->BBB)
    {
      rv = Py_True;
    }
    else
    {
      rv = Py_False;
    }
    break;
  case Py_NE:
    if (self->BBB != pyother->BBB)
    {
      rv = Py_True;
    }
    else
    {
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
PyObject* rv = Py_NotImplemented;
if (PyObject_IsInstance(other, (PyObject*) &PY_DataBuffer_Type))
{
  PY_DataBuffer* pyother = (PY_DataBuffer*) other;
  switch (opid)
  {
  case Py_EQ:
    if (self->BBB == pyother->BBB)
    {
      rv = Py_True;
    }
    else
    {
      rv = Py_False;
    }
    break;
  case Py_NE:
    if (self->BBB != pyother->BBB)
    {
      rv = Py_True;
    }
    else
    {
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
PyObject* rv = Py_NotImplemented;
if (PyObject_IsInstance(other, (PyObject*) &PY_DataView_Type))
{
  PY_DataView* pyother = (PY_DataView*) other;
  switch (opid)
  {
  case Py_EQ:
    if (self->BBB == pyother->BBB)
    {
      rv = Py_True;
    }
    else
    {
      rv = Py_False;
    }
    break;
  case Py_NE:
    if (self->BBB != pyother->BBB)
    {
      rv = Py_True;
    }
    else
    {
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


// ----------------------------------------------------------------------

// splicer begin function.name_is_valid
// Allow any sort of object.  This will allow None to return False.
PyObject* name;
const char* kwcpp = "name";
char* kw_list[] = { (char*) kwcpp+0, NULL };
bool rv = false;

if (!PyArg_ParseTupleAndKeywords(args, kwds, "O:nameIsValid", kw_list,
                                 &name))
{
  return NULL;
}
if (PyString_Check(name))
{
  rv = nameIsValid(PyString_AS_STRING(name));
}
return PyBool_FromLong(rv);
// splicer end function.name_is_valid
