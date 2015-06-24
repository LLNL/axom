// pyDataStoretype.cpp
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
// splicer begin class.DataStore.include
// splicer end class.DataStore.include
namespace asctoolkit {
namespace sidre {
// splicer begin class.DataStore.C_definition
// splicer end class.DataStore.C_definition
// splicer begin class.DataStore.additional_methods
// splicer end class.DataStore.additional_methods
static int
PY_DataStore_tp_init (PY_DataStore *self, PyObject *args, PyObject *kwds)
{
// splicer begin class.DataStore.type.init
DataStore * ds = new DataStore();
self->BBB = ds;
return 0;
// splicer end class.DataStore.type.init
}

static char PY_datastore_get_root__doc__[] =
"documentation"
;

static PyObject *
PY_datastore_get_root(
  PY_DataStore *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataStore.method.getRoot
    DataGroup * rv = self->BBB->getRoot();
    return Py_BuildValue("O&", PP_DataGroup_to_Object, rv);
// splicer end class.DataStore.method.getRoot
}

static char PY_datastore_get_buffer__doc__[] =
"documentation"
;

static PyObject *
PY_datastore_get_buffer(
  PY_DataStore *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataStore.method.getBuffer
    ATK_IndexType idx;
    const char *kwcpp = "idx";
    char *kw_list[] = { (char *) kwcpp+0 };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O:getBuffer", kw_list,
        &idx))
    {
        return NULL;
    }
    DataBuffer * rv = self->BBB->getBuffer(idx);
    return Py_BuildValue("O&", PP_DataBuffer_to_Object, rv);
// splicer end class.DataStore.method.getBuffer
}

static char PY_datastore_create_buffer__doc__[] =
"documentation"
;

static PyObject *
PY_datastore_create_buffer(
  PY_DataStore *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataStore.method.createBuffer
    DataBuffer * rv = self->BBB->createBuffer();
    return Py_BuildValue("O&", PP_DataBuffer_to_Object, rv);
// splicer end class.DataStore.method.createBuffer
}

static char PY_datastore_destroy_buffer__doc__[] =
"documentation"
;

static PyObject *
PY_datastore_destroy_buffer(
  PY_DataStore *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataStore.method.destroyBuffer
    ATK_IndexType id;
    const char *kwcpp = "id";
    char *kw_list[] = { (char *) kwcpp+0 };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O:destroyBuffer", kw_list,
        &id))
    {
        return NULL;
    }
    self->BBB->destroyBuffer(id);
    Py_RETURN_NONE;
// splicer end class.DataStore.method.destroyBuffer
}

static char PY_datastore_get_num_buffers__doc__[] =
"documentation"
;

static PyObject *
PY_datastore_get_num_buffers(
  PY_DataStore *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataStore.method.getNumBuffers
    size_t rv = self->BBB->getNumBuffers();
    return Py_BuildValue("O", &rv);
// splicer end class.DataStore.method.getNumBuffers
}

static char PY_datastore_print__doc__[] =
"documentation"
;

static PyObject *
PY_datastore_print(
  PY_DataStore *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataStore.method.print
    self->BBB->print();
    Py_RETURN_NONE;
// splicer end class.DataStore.method.print
}
static PyMethodDef PY_DataStore_methods[] = {
{"getRoot", (PyCFunction)PY_datastore_get_root, METH_NOARGS, PY_datastore_get_root__doc__},
{"getBuffer", (PyCFunction)PY_datastore_get_buffer, METH_VARARGS|METH_KEYWORDS, PY_datastore_get_buffer__doc__},
{"createBuffer", (PyCFunction)PY_datastore_create_buffer, METH_NOARGS, PY_datastore_create_buffer__doc__},
{"destroyBuffer", (PyCFunction)PY_datastore_destroy_buffer, METH_VARARGS|METH_KEYWORDS, PY_datastore_destroy_buffer__doc__},
{"getNumBuffers", (PyCFunction)PY_datastore_get_num_buffers, METH_NOARGS, PY_datastore_get_num_buffers__doc__},
{"print", (PyCFunction)PY_datastore_print, METH_NOARGS, PY_datastore_print__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

static char DataStore__doc__[] =
"virtual class"
;

/* static */
PyTypeObject PY_DataStore_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "sidre.DataStore",                       /* tp_name */
        sizeof(PY_DataStore),         /* tp_basicsize */
        0,                              /* tp_itemsize */
        /* Methods to implement standard operations */
        (destructor)0,                 /* tp_dealloc */
        (printfunc)0,                   /* tp_print */
        (getattrfunc)0,                 /* tp_getattr */
        (setattrfunc)0,                 /* tp_setattr */
#ifdef IS_PY3K
        0,                               /* tp_reserved */
#else
        (cmpfunc)0,                     /* tp_compare */
#endif
        (reprfunc)0,                    /* tp_repr */
        /* Method suites for standard classes */
        0,                              /* tp_as_number */
        0,                              /* tp_as_sequence */
        0,                              /* tp_as_mapping */
        /* More standard operations (here for binary compatibility) */
        (hashfunc)0,                    /* tp_hash */
        (ternaryfunc)0,                 /* tp_call */
        (reprfunc)0,                    /* tp_str */
        (getattrofunc)0,                /* tp_getattro */
        (setattrofunc)0,                /* tp_setattro */
        /* Functions to access object as input/output buffer */
        0,                              /* tp_as_buffer */
        /* Flags to define presence of optional/expanded features */
        Py_TPFLAGS_DEFAULT,             /* tp_flags */
        DataStore__doc__,         /* tp_doc */
        /* Assigned meaning in release 2.0 */
        /* call function for all accessible objects */
        (traverseproc)0,                /* tp_traverse */
        /* delete references to contained objects */
        (inquiry)0,                     /* tp_clear */
        /* Assigned meaning in release 2.1 */
        /* rich comparisons */
        (richcmpfunc)0,                 /* tp_richcompare */
        /* weak reference enabler */
        0,                              /* tp_weaklistoffset */
        /* Added in release 2.2 */
        /* Iterators */
        (getiterfunc)0,                 /* tp_iter */
        (iternextfunc)0,                /* tp_iternext */
        /* Attribute descriptor and subclassing stuff */
        PY_DataStore_methods,                             /* tp_methods */
        0,                              /* tp_members */
        0,                             /* tp_getset */
        0,                              /* tp_base */
        0,                              /* tp_dict */
        (descrgetfunc)0,                /* tp_descr_get */
        (descrsetfunc)0,                /* tp_descr_set */
        0,                              /* tp_dictoffset */
        (initproc)PY_DataStore_tp_init,                   /* tp_init */
        (allocfunc)0,                  /* tp_alloc */
        (newfunc)0,                    /* tp_new */
        (freefunc)0,                   /* tp_free */
        (inquiry)0,                     /* tp_is_gc */
        0,                              /* tp_bases */
        0,                              /* tp_mro */
        0,                              /* tp_cache */
        0,                              /* tp_subclasses */
        0,                              /* tp_weaklist */
        (destructor)0,                 /* tp_del */
        0,                              /* tp_version_tag */
#ifdef IS_PY3K
        (destructor)0,                  /* tp_finalize */
#endif
};

}  // namespace asctoolkit
}  // namespace sidre
