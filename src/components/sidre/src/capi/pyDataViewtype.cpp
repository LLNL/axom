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
// splicer begin class.DataView.impl.include
// splicer end class.DataView.impl.include
namespace asctoolkit {
namespace sidre {
// splicer begin class.DataView.impl.C_definition
// splicer end class.DataView.impl.C_definition
// splicer begin class.DataView.impl.additional_methods
// splicer end class.DataView.impl.additional_methods
static int
PY_DataView_tp_init (PY_DataView *self, PyObject *args, PyObject *kwds)
{
// splicer begin class.DataView.type.init
    PyErr_SetString(PyExc_NotImplementedError, "init");
    return -1;
// splicer end class.DataView.type.init
}
static PyObject *
PY_DataView_tp_richcompare (PY_DataView *self, PyObject *other, int opid)
{
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
}

static char PY_dataview_allocate_simple__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_allocate_simple(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.allocate_simple
    self->BBB->allocate();
    Py_RETURN_NONE;
// splicer end class.DataView.method.allocate_simple
}

static char PY_dataview_allocate_from_type__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_allocate_from_type(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.allocate_from_type
    int type;
    SIDRE_SidreLength num_elems;
    const char *kwcpp = "type\0num_elems";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "il:allocate", kw_list,
        &type, &num_elems))
    {
        return NULL;
    }
    self->BBB->allocate(getTypeID(type), num_elems);
    Py_RETURN_NONE;
// splicer end class.DataView.method.allocate_from_type
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
    SIDRE_SidreLength num_elems;
    const char *kwcpp = "num_elems";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "l:reallocate", kw_list,
        &num_elems))
    {
        return NULL;
    }
    self->BBB->reallocate(num_elems);
    Py_RETURN_NONE;
// splicer end class.DataView.method.reallocate
}

static char PY_dataview_apply_0__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_apply_0(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.apply_0
    self->BBB->apply();
    Py_RETURN_NONE;
// splicer end class.DataView.method.apply_0
}

static char PY_dataview_attach_buffer__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_attach_buffer(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.attach_buffer
    PY_DataBuffer * buff;
    DataBuffer * buff_ptr;
    const char *kwcpp = "buff";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!:attachBuffer", kw_list,
        &PY_DataBuffer_Type, &buff))
    {
        return NULL;
    }
    buff_ptr = (buff ? buff->BBB : NULL);
    self->BBB->attachBuffer(buff_ptr);
    Py_RETURN_NONE;
// splicer end class.DataView.method.attach_buffer
}

static char PY_dataview_apply_nelems__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_apply_nelems(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.apply_nelems
    SIDRE_SidreLength num_elems;
    const char *kwcpp = "num_elems";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "l:apply", kw_list,
        &num_elems))
    {
        return NULL;
    }
    self->BBB->apply(num_elems);
    Py_RETURN_NONE;
// splicer end class.DataView.method.apply_nelems
}

static char PY_dataview_apply_nelems_offset__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_apply_nelems_offset(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.apply_nelems_offset
    SIDRE_SidreLength num_elems;
    SIDRE_SidreLength offset;
    const char *kwcpp = "num_elems\0offset";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+10, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll:apply", kw_list,
        &num_elems, &offset))
    {
        return NULL;
    }
    self->BBB->apply(num_elems, offset);
    Py_RETURN_NONE;
// splicer end class.DataView.method.apply_nelems_offset
}

static char PY_dataview_apply_nelems_offset_stride__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_apply_nelems_offset_stride(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.apply_nelems_offset_stride
    SIDRE_SidreLength num_elems;
    SIDRE_SidreLength offset;
    SIDRE_SidreLength stride;
    const char *kwcpp = "num_elems\0offset\0stride";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+10,(char *) kwcpp+17, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "lll:apply", kw_list,
        &num_elems, &offset, &stride))
    {
        return NULL;
    }
    self->BBB->apply(num_elems, offset, stride);
    Py_RETURN_NONE;
// splicer end class.DataView.method.apply_nelems_offset_stride
}

static char PY_dataview_apply_type_nelems__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_apply_type_nelems(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.apply_type_nelems
    int type;
    SIDRE_SidreLength num_elems;
    const char *kwcpp = "type\0num_elems";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "il:apply", kw_list,
        &type, &num_elems))
    {
        return NULL;
    }
    self->BBB->apply(getTypeID(type), num_elems);
    Py_RETURN_NONE;
// splicer end class.DataView.method.apply_type_nelems
}

static char PY_dataview_apply_type_nelems_offset__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_apply_type_nelems_offset(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.apply_type_nelems_offset
    int type;
    SIDRE_SidreLength num_elems;
    SIDRE_SidreLength offset;
    const char *kwcpp = "type\0num_elems\0offset";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5,(char *) kwcpp+15, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ill:apply", kw_list,
        &type, &num_elems, &offset))
    {
        return NULL;
    }
    self->BBB->apply(getTypeID(type), num_elems, offset);
    Py_RETURN_NONE;
// splicer end class.DataView.method.apply_type_nelems_offset
}

static char PY_dataview_apply_type_nelems_offset_stride__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_apply_type_nelems_offset_stride(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.apply_type_nelems_offset_stride
    int type;
    SIDRE_SidreLength num_elems;
    SIDRE_SidreLength offset;
    SIDRE_SidreLength stride;
    const char *kwcpp = "type\0num_elems\0offset\0stride";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5,(char *) kwcpp+15,(char *) kwcpp+22, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "illl:apply", kw_list,
        &type, &num_elems, &offset, &stride))
    {
        return NULL;
    }
    self->BBB->apply(getTypeID(type), num_elems, offset, stride);
    Py_RETURN_NONE;
// splicer end class.DataView.method.apply_type_nelems_offset_stride
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
// splicer begin class.DataView.method.has_buffer
    bool rv = self->BBB->hasBuffer();
    return PyBool_FromLong(rv);
// splicer end class.DataView.method.has_buffer
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
// splicer begin class.DataView.method.is_opaque
    bool rv = self->BBB->isOpaque();
    return PyBool_FromLong(rv);
// splicer end class.DataView.method.is_opaque
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
// splicer begin class.DataView.method.get_name
    const std::string & rv = self->BBB->getName();
    return PyString_FromString(rv.c_str());
// splicer end class.DataView.method.get_name
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
// splicer begin class.DataView.method.get_opaque
    void * rv = self->BBB->getOpaque();
    return Py_BuildValue("O", rv);
// splicer end class.DataView.method.get_opaque
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
// splicer begin class.DataView.method.get_buffer
    DataBuffer * rv = self->BBB->getBuffer();
    PY_DataBuffer * rv_obj = PyObject_New(PY_DataBuffer, &PY_DataBuffer_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataView.method.get_buffer
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
// splicer begin class.DataView.method.get_data_pointer
    void * rv = self->BBB->getDataPointer();
    return Py_BuildValue("O", rv);
// splicer end class.DataView.method.get_data_pointer
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
// splicer begin class.DataView.method.get_owning_group
    DataGroup * rv = self->BBB->getOwningGroup();
    PY_DataGroup * rv_obj = PyObject_New(PY_DataGroup, &PY_DataGroup_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataView.method.get_owning_group
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
// splicer begin class.DataView.method.get_type_id
    TypeID rv = self->BBB->getTypeID();
    return Py_BuildValue("i", rv);
// splicer end class.DataView.method.get_type_id
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
// splicer begin class.DataView.method.get_total_bytes
    size_t rv = self->BBB->getTotalBytes();
    return PyInt_FromLong(rv);
// splicer end class.DataView.method.get_total_bytes
}

static char PY_dataview_get_num_elements__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_get_num_elements(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.method.get_num_elements
    size_t rv = self->BBB->getNumElements();
    return PyInt_FromLong(rv);
// splicer end class.DataView.method.get_num_elements
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

static char PY_dataview_allocate__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_allocate(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.allocate
    int numNamedArgs = (kwds ? PyDict_Size(kwds) : 0);
    int numArgs = PyTuple_GET_SIZE(args);
    int totArgs = numArgs + numNamedArgs;
    PyObject *rvobj;
    {
        rvobj = PY_dataview_allocate_simple(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    {
        rvobj = PY_dataview_allocate_from_type(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
    return NULL;
// splicer end class.DataView.allocate
}

static char PY_dataview_apply__doc__[] =
"documentation"
;

static PyObject *
PY_dataview_apply(
  PY_DataView *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataView.apply
    int numNamedArgs = (kwds ? PyDict_Size(kwds) : 0);
    int numArgs = PyTuple_GET_SIZE(args);
    int totArgs = numArgs + numNamedArgs;
    PyObject *rvobj;
    {
        rvobj = PY_dataview_apply_0(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    {
        rvobj = PY_dataview_apply_nelems(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    {
        rvobj = PY_dataview_apply_nelems_offset(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    {
        rvobj = PY_dataview_apply_nelems_offset_stride(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    {
        rvobj = PY_dataview_apply_type_nelems(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    {
        rvobj = PY_dataview_apply_type_nelems_offset(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    {
        rvobj = PY_dataview_apply_type_nelems_offset_stride(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
    return NULL;
// splicer end class.DataView.apply
}
// splicer begin class.DataView.impl.after_methods
// splicer end class.DataView.impl.after_methods
static PyMethodDef PY_DataView_methods[] = {
{"allocate_simple", (PyCFunction)PY_dataview_allocate_simple, METH_NOARGS, PY_dataview_allocate_simple__doc__},
{"allocate_from_type", (PyCFunction)PY_dataview_allocate_from_type, METH_VARARGS|METH_KEYWORDS, PY_dataview_allocate_from_type__doc__},
{"reallocate", (PyCFunction)PY_dataview_reallocate, METH_VARARGS|METH_KEYWORDS, PY_dataview_reallocate__doc__},
{"apply_0", (PyCFunction)PY_dataview_apply_0, METH_NOARGS, PY_dataview_apply_0__doc__},
{"attachBuffer", (PyCFunction)PY_dataview_attach_buffer, METH_VARARGS|METH_KEYWORDS, PY_dataview_attach_buffer__doc__},
{"apply_nelems", (PyCFunction)PY_dataview_apply_nelems, METH_VARARGS|METH_KEYWORDS, PY_dataview_apply_nelems__doc__},
{"apply_nelems_offset", (PyCFunction)PY_dataview_apply_nelems_offset, METH_VARARGS|METH_KEYWORDS, PY_dataview_apply_nelems_offset__doc__},
{"apply_nelems_offset_stride", (PyCFunction)PY_dataview_apply_nelems_offset_stride, METH_VARARGS|METH_KEYWORDS, PY_dataview_apply_nelems_offset_stride__doc__},
{"apply_type_nelems", (PyCFunction)PY_dataview_apply_type_nelems, METH_VARARGS|METH_KEYWORDS, PY_dataview_apply_type_nelems__doc__},
{"apply_type_nelems_offset", (PyCFunction)PY_dataview_apply_type_nelems_offset, METH_VARARGS|METH_KEYWORDS, PY_dataview_apply_type_nelems_offset__doc__},
{"apply_type_nelems_offset_stride", (PyCFunction)PY_dataview_apply_type_nelems_offset_stride, METH_VARARGS|METH_KEYWORDS, PY_dataview_apply_type_nelems_offset_stride__doc__},
{"hasBuffer", (PyCFunction)PY_dataview_has_buffer, METH_NOARGS, PY_dataview_has_buffer__doc__},
{"isOpaque", (PyCFunction)PY_dataview_is_opaque, METH_NOARGS, PY_dataview_is_opaque__doc__},
{"getName", (PyCFunction)PY_dataview_get_name, METH_NOARGS, PY_dataview_get_name__doc__},
{"getOpaque", (PyCFunction)PY_dataview_get_opaque, METH_NOARGS, PY_dataview_get_opaque__doc__},
{"getBuffer", (PyCFunction)PY_dataview_get_buffer, METH_NOARGS, PY_dataview_get_buffer__doc__},
{"getDataPointer", (PyCFunction)PY_dataview_get_data_pointer, METH_NOARGS, PY_dataview_get_data_pointer__doc__},
{"getOwningGroup", (PyCFunction)PY_dataview_get_owning_group, METH_NOARGS, PY_dataview_get_owning_group__doc__},
{"getTypeID", (PyCFunction)PY_dataview_get_type_id, METH_NOARGS, PY_dataview_get_type_id__doc__},
{"getTotalBytes", (PyCFunction)PY_dataview_get_total_bytes, METH_NOARGS, PY_dataview_get_total_bytes__doc__},
{"getNumElements", (PyCFunction)PY_dataview_get_num_elements, METH_NOARGS, PY_dataview_get_num_elements__doc__},
{"print", (PyCFunction)PY_dataview_print, METH_NOARGS, PY_dataview_print__doc__},
{"allocate", (PyCFunction)PY_dataview_allocate, METH_VARARGS|METH_KEYWORDS, PY_dataview_allocate__doc__},
{"apply", (PyCFunction)PY_dataview_apply, METH_VARARGS|METH_KEYWORDS, PY_dataview_apply__doc__},
// splicer begin class.DataView.PyMethodDef
// splicer end class.DataView.PyMethodDef
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
        (richcmpfunc)PY_DataView_tp_richcompare,                 /* tp_richcompare */
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
