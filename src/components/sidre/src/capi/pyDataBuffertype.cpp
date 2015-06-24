// pyDataBuffertype.cpp
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
// splicer begin class.DataBuffer.include
// splicer end class.DataBuffer.include
namespace asctoolkit {
namespace sidre {
// splicer begin class.DataBuffer.C_definition
// splicer end class.DataBuffer.C_definition
// splicer begin class.DataBuffer.additional_methods
// splicer end class.DataBuffer.additional_methods
static int
PY_DataBuffer_tp_init (PY_DataBuffer *self, PyObject *args, PyObject *kwds)
{
// splicer begin class.DataBuffer.type.init
    PyErr_SetString(PyExc_NotImplementedError, "init");
    return -1;
// splicer end class.DataBuffer.type.init
}

static char PY_databuffer_get_index__doc__[] =
"documentation"
;

static PyObject *
PY_databuffer_get_index(
  PY_DataBuffer *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataBuffer.method.getIndex
    IndexType rv = self->BBB->getIndex();
    return Py_BuildValue("O", &rv);
// splicer end class.DataBuffer.method.getIndex
}

static char PY_databuffer_get_num_views__doc__[] =
"documentation"
;

static PyObject *
PY_databuffer_get_num_views(
  PY_DataBuffer *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataBuffer.method.getNumViews
    size_t rv = self->BBB->getNumViews();
    return Py_BuildValue("O", &rv);
// splicer end class.DataBuffer.method.getNumViews
}

static char PY_databuffer_declare__doc__[] =
"documentation"
;

static PyObject *
PY_databuffer_declare(
  PY_DataBuffer *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataBuffer.method.declare
    int type;
    ATK_SidreLength len;
    const char *kwcpp = "type\0len";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5 };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO:declare", kw_list,
        &type, &len))
    {
        return NULL;
    }
    self->BBB->declare(type, len);
    Py_RETURN_NONE;
// splicer end class.DataBuffer.method.declare
}

static char PY_databuffer_declare_external__doc__[] =
"documentation"
;

static PyObject *
PY_databuffer_declare_external(
  PY_DataBuffer *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataBuffer.method.declareExternal
    void * external_data;
    int type;
    ATK_SidreLength len;
    const char *kwcpp = "external_data\0type\0len";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+14,(char *) kwcpp+19 };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOO:declareExternal", kw_list,
        &external_data, &type, &len))
    {
        return NULL;
    }
    self->BBB->declareExternal(external_data, type, len);
    Py_RETURN_NONE;
// splicer end class.DataBuffer.method.declareExternal
}

static char PY_databuffer_allocate_existing__doc__[] =
"documentation"
;

static PyObject *
PY_databuffer_allocate_existing(
  PY_DataBuffer *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataBuffer.method.allocate
    self->BBB->allocate();
    Py_RETURN_NONE;
// splicer end class.DataBuffer.method.allocate
}

static char PY_databuffer_allocate_from_type__doc__[] =
"documentation"
;

static PyObject *
PY_databuffer_allocate_from_type(
  PY_DataBuffer *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataBuffer.method.allocate
    int type;
    ATK_SidreLength len;
    const char *kwcpp = "type\0len";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5 };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO:allocate", kw_list,
        &type, &len))
    {
        return NULL;
    }
    self->BBB->allocate(type, len);
    Py_RETURN_NONE;
// splicer end class.DataBuffer.method.allocate
}

static char PY_databuffer_reallocate__doc__[] =
"documentation"
;

static PyObject *
PY_databuffer_reallocate(
  PY_DataBuffer *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataBuffer.method.reallocate
    int type;
    ATK_SidreLength len;
    const char *kwcpp = "type\0len";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5 };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO:reallocate", kw_list,
        &type, &len))
    {
        return NULL;
    }
    self->BBB->reallocate(type, len);
    Py_RETURN_NONE;
// splicer end class.DataBuffer.method.reallocate
}

static char PY_databuffer_is_external__doc__[] =
"documentation"
;

static PyObject *
PY_databuffer_is_external(
  PY_DataBuffer *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataBuffer.method.isExternal
    bool rv = self->BBB->isExternal();
    return Py_BuildValue("O", &rv);
// splicer end class.DataBuffer.method.isExternal
}

static char PY_databuffer_get_data__doc__[] =
"documentation"
;

static PyObject *
PY_databuffer_get_data(
  PY_DataBuffer *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataBuffer.method.getData
    void * rv = self->BBB->getData();
    return Py_BuildValue("O", rv);
// splicer end class.DataBuffer.method.getData
}

static char PY_databuffer_get_total_bytes__doc__[] =
"documentation"
;

static PyObject *
PY_databuffer_get_total_bytes(
  PY_DataBuffer *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataBuffer.method.getTotalBytes
    size_t rv = self->BBB->getTotalBytes();
    return Py_BuildValue("O", &rv);
// splicer end class.DataBuffer.method.getTotalBytes
}
static PyMethodDef PY_DataBuffer_methods[] = {
{"getIndex", (PyCFunction)PY_databuffer_get_index, METH_NOARGS, PY_databuffer_get_index__doc__},
{"getNumViews", (PyCFunction)PY_databuffer_get_num_views, METH_NOARGS, PY_databuffer_get_num_views__doc__},
{"declare", (PyCFunction)PY_databuffer_declare, METH_VARARGS|METH_KEYWORDS, PY_databuffer_declare__doc__},
{"declareExternal", (PyCFunction)PY_databuffer_declare_external, METH_VARARGS|METH_KEYWORDS, PY_databuffer_declare_external__doc__},
{"allocate_existing", (PyCFunction)PY_databuffer_allocate_existing, METH_NOARGS, PY_databuffer_allocate_existing__doc__},
{"allocate_from_type", (PyCFunction)PY_databuffer_allocate_from_type, METH_VARARGS|METH_KEYWORDS, PY_databuffer_allocate_from_type__doc__},
{"reallocate", (PyCFunction)PY_databuffer_reallocate, METH_VARARGS|METH_KEYWORDS, PY_databuffer_reallocate__doc__},
{"isExternal", (PyCFunction)PY_databuffer_is_external, METH_NOARGS, PY_databuffer_is_external__doc__},
{"getData", (PyCFunction)PY_databuffer_get_data, METH_NOARGS, PY_databuffer_get_data__doc__},
{"getTotalBytes", (PyCFunction)PY_databuffer_get_total_bytes, METH_NOARGS, PY_databuffer_get_total_bytes__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

static char DataBuffer__doc__[] =
"virtual class"
;

/* static */
PyTypeObject PY_DataBuffer_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "sidre.DataBuffer",                       /* tp_name */
        sizeof(PY_DataBuffer),         /* tp_basicsize */
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
        DataBuffer__doc__,         /* tp_doc */
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
        PY_DataBuffer_methods,                             /* tp_methods */
        0,                              /* tp_members */
        0,                             /* tp_getset */
        0,                              /* tp_base */
        0,                              /* tp_dict */
        (descrgetfunc)0,                /* tp_descr_get */
        (descrsetfunc)0,                /* tp_descr_set */
        0,                              /* tp_dictoffset */
        (initproc)PY_DataBuffer_tp_init,                   /* tp_init */
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
