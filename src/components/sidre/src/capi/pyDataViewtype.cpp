// pyDataViewtype.cpp
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
// splicer begin class.DataView.include
// splicer end class.DataView.include
namespace asctoolkit {
namespace sidre {
// splicer begin class.DataView.C_definition
// splicer end class.DataView.C_definition
// splicer begin class.DataView.additional_methods
// splicer end class.DataView.additional_methods
static int
PY_DataView_tp_init (PY_DataView *self, PyObject *args, PyObject *kwds)
{
// splicer begin class.DataView.type.init
    PyErr_SetString(PyExc_NotImplementedError, "init");
    return -1;
// splicer end class.DataView.type.init
}

static char PY_dataview_declare__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_declare(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.declare
    int type;
    ATK_SidreLength len;
    const char *kwcpp = "type\0len";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "il:declare", kw_list,
        &type, &len))
    {
        return NULL;
    }
    self->BBB->declare(getTypeID(type), len);
    Py_RETURN_NONE;
// splicer end class.DataView.method.declare
}

static char PY_dataview_allocate__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_allocate(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.allocate
    int type;
    ATK_SidreLength len;
    const char *kwcpp = "type\0len";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "il:allocate", kw_list,
        &type, &len))
    {
        return NULL;
    }
    self->BBB->allocate(getTypeID(type), len);
    Py_RETURN_NONE;
// splicer end class.DataView.method.allocate
}

static char PY_dataview_reallocate__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_reallocate(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.reallocate
    int type;
    ATK_SidreLength len;
    const char *kwcpp = "type\0len";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "il:reallocate", kw_list,
        &type, &len))
    {
        return NULL;
    }
    self->BBB->reallocate(getTypeID(type), len);
    Py_RETURN_NONE;
// splicer end class.DataView.method.reallocate
}

static char PY_dataview_has_buffer__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_has_buffer(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.hasBuffer
    bool rv = self->BBB->hasBuffer();
    return Py_BuildValue("O", &rv);
// splicer end class.DataView.method.hasBuffer
}

static char PY_dataview_is_opaque__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_is_opaque(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.isOpaque
    bool rv = self->BBB->isOpaque();
    return Py_BuildValue("O", &rv);
// splicer end class.DataView.method.isOpaque
}

static char PY_dataview_get_name__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_get_name(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.getName
    const std::string & rv = self->BBB->getName();
    return Py_BuildValue("s", rv.c_str());
// splicer end class.DataView.method.getName
}

static char PY_dataview_get_opaque__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_get_opaque(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.getOpaque
    void * rv = self->BBB->getOpaque();
    return Py_BuildValue("O", rv);
// splicer end class.DataView.method.getOpaque
}

static char PY_dataview_get_buffer__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_get_buffer(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.getBuffer
    DataBuffer * rv = self->BBB->getBuffer();
    PY_DataBuffer * rv_obj = PyObject_New(PY_DataBuffer, &PY_DataBuffer_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataView.method.getBuffer
}

static char PY_dataview_get_data_pointer__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_get_data_pointer(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.getDataPointer
    void * rv = self->BBB->getDataPointer();
    return Py_BuildValue("O", rv);
// splicer end class.DataView.method.getDataPointer
}

static char PY_dataview_get_owning_group__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_get_owning_group(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.getOwningGroup
    DataGroup * rv = self->BBB->getOwningGroup();
    PY_DataGroup * rv_obj = PyObject_New(PY_DataGroup, &PY_DataGroup_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataView.method.getOwningGroup
}

static char PY_dataview_get_type_id__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_get_type_id(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.getTypeID
    TypeID rv = self->BBB->getTypeID();
    return Py_BuildValue("i", &rv);
// splicer end class.DataView.method.getTypeID
}

static char PY_dataview_get_total_bytes__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_get_total_bytes(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.getTotalBytes
    size_t rv = self->BBB->getTotalBytes();
    return Py_BuildValue("O", &rv);
// splicer end class.DataView.method.getTotalBytes
}

static char PY_dataview_get_number_of_elements__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_get_number_of_elements(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.getNumberOfElements
    size_t rv = self->BBB->getNumberOfElements();
    return Py_BuildValue("O", &rv);
// splicer end class.DataView.method.getNumberOfElements
}

static char PY_dataview_print__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_print(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.print
    self->BBB->print();
    Py_RETURN_NONE;
// splicer end class.DataView.method.print
}
static PyMethodDef PY_DataView_methods[] = {
{"declare", (PyCFunction)PY_dataview_declare, METH_VARARGS|METH_KEYWORDS, PY_dataview_declare__doc__},
{"allocate", (PyCFunction)PY_dataview_allocate, METH_VARARGS|METH_KEYWORDS, PY_dataview_allocate__doc__},
{"reallocate", (PyCFunction)PY_dataview_reallocate, METH_VARARGS|METH_KEYWORDS, PY_dataview_reallocate__doc__},
{"hasBuffer", (PyCFunction)PY_dataview_has_buffer, METH_NOARGS, PY_dataview_has_buffer__doc__},
{"isOpaque", (PyCFunction)PY_dataview_is_opaque, METH_NOARGS, PY_dataview_is_opaque__doc__},
{"getName", (PyCFunction)PY_dataview_get_name, METH_NOARGS, PY_dataview_get_name__doc__},
{"getOpaque", (PyCFunction)PY_dataview_get_opaque, METH_NOARGS, PY_dataview_get_opaque__doc__},
{"getBuffer", (PyCFunction)PY_dataview_get_buffer, METH_NOARGS, PY_dataview_get_buffer__doc__},
{"getDataPointer", (PyCFunction)PY_dataview_get_data_pointer, METH_NOARGS, PY_dataview_get_data_pointer__doc__},
{"getOwningGroup", (PyCFunction)PY_dataview_get_owning_group, METH_NOARGS, PY_dataview_get_owning_group__doc__},
{"getTypeID", (PyCFunction)PY_dataview_get_type_id, METH_NOARGS, PY_dataview_get_type_id__doc__},
{"getTotalBytes", (PyCFunction)PY_dataview_get_total_bytes, METH_NOARGS, PY_dataview_get_total_bytes__doc__},
{"getNumberOfElements", (PyCFunction)PY_dataview_get_number_of_elements, METH_NOARGS, PY_dataview_get_number_of_elements__doc__},
{"print", (PyCFunction)PY_dataview_print, METH_NOARGS, PY_dataview_print__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

static char DataView__doc__[] =
"virtual class"
;

/* static */
PyTypeObject PY_DataView_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "sidre.DataView",                       /* tp_name */
        sizeof(PY_DataView),         /* tp_basicsize */
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
        DataView__doc__,         /* tp_doc */
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
        PY_DataView_methods,                             /* tp_methods */
        0,                              /* tp_members */
        0,                             /* tp_getset */
        0,                              /* tp_base */
        0,                              /* tp_dict */
        (descrgetfunc)0,                /* tp_descr_get */
        (descrsetfunc)0,                /* tp_descr_set */
        0,                              /* tp_dictoffset */
        (initproc)PY_DataView_tp_init,                   /* tp_init */
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
