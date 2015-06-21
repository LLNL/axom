// pyExClass1type.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
#include "pyUserLibrarymodule.hpp"
// splicer begin class.ExClass1.include
// splicer end class.ExClass1.include
namespace example {
namespace nested {
// splicer begin class.ExClass1.C_definition
// splicer end class.ExClass1.C_definition
// splicer begin class.ExClass1.additional_methods
// splicer end class.ExClass1.additional_methods
static PyObject *
PP_ExClass1_tp_repr (PP_ExClass1 *self)
{
// splicer begin class.ExClass1.type.repr
    repr code
// splicer end class.ExClass1.type.repr
}
static int
PP_ExClass1_tp_init (PP_ExClass1 *self, PyObject *args, PyObject *kwds)
{
// splicer begin class.ExClass1.type.init
    init code
// splicer end class.ExClass1.type.init
}

static char PP_exclass1_new__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_new(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.new
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass1.method.new
}

static char PP_exclass1_delete__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_delete(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.delete
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass1.method.delete
}

static char PP_exclass1_increment_count__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_increment_count(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.incrementCount
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass1.method.incrementCount
}

static char PP_exclass1_get_name__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_get_name(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.getName
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass1.method.getName
}

static char PP_exclass1_get_name_length__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_get_name_length(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.get_name_length
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass1.method.get_name_length
}

static char PP_exclass1_get_root__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_get_root(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.getRoot
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass1.method.getRoot
}

static char PP_exclass1_get_value_from_int__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_get_value_from_int(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.getValue
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass1.method.getValue
}

static char PP_exclass1_get_value_1__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_get_value_1(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.getValue
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass1.method.getValue
}

static char PP_exclass1_get_addr__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_get_addr(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.getAddr
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass1.method.getAddr
}

static char PP_exclass1_has_addr__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_has_addr(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.hasAddr
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass1.method.hasAddr
}

static char PP_exclass1_splicer_special__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_splicer_special(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.SplicerSpecial
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass1.method.SplicerSpecial
}
static PyMethodDef PP_ExClass1_methods[] = {
{"new", (PyCFunction)PP_exclass1_new, METH_VARARGS|METH_KEYWORDS, PP_exclass1_new__doc__},
{"delete", (PyCFunction)PP_exclass1_delete, METH_NOARGS, PP_exclass1_delete__doc__},
{"incrementCount", (PyCFunction)PP_exclass1_increment_count, METH_VARARGS|METH_KEYWORDS, PP_exclass1_increment_count__doc__},
{"getName", (PyCFunction)PP_exclass1_get_name, METH_NOARGS, PP_exclass1_get_name__doc__},
{"get_name_length", (PyCFunction)PP_exclass1_get_name_length, METH_NOARGS, PP_exclass1_get_name_length__doc__},
{"getRoot", (PyCFunction)PP_exclass1_get_root, METH_NOARGS, PP_exclass1_get_root__doc__},
{"getValue_from_int", (PyCFunction)PP_exclass1_get_value_from_int, METH_VARARGS|METH_KEYWORDS, PP_exclass1_get_value_from_int__doc__},
{"getValue_1", (PyCFunction)PP_exclass1_get_value_1, METH_VARARGS|METH_KEYWORDS, PP_exclass1_get_value_1__doc__},
{"getAddr", (PyCFunction)PP_exclass1_get_addr, METH_NOARGS, PP_exclass1_get_addr__doc__},
{"hasAddr", (PyCFunction)PP_exclass1_has_addr, METH_VARARGS|METH_KEYWORDS, PP_exclass1_has_addr__doc__},
{"SplicerSpecial", (PyCFunction)PP_exclass1_splicer_special, METH_NOARGS, PP_exclass1_splicer_special__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

static char ExClass1__doc__[] =
"virtual class"
;

/* static */
PyTypeObject PP_ExClass1_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "userlibrary.ExClass1",                       /* tp_name */
        sizeof(PP_ExClass1),         /* tp_basicsize */
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
        (reprfunc)PP_ExClass1_tp_repr,                    /* tp_repr */
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
        ExClass1__doc__,         /* tp_doc */
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
        PP_ExClass1_methods,                             /* tp_methods */
        0,                              /* tp_members */
        0,                             /* tp_getset */
        0,                              /* tp_base */
        0,                              /* tp_dict */
        (descrgetfunc)0,                /* tp_descr_get */
        (descrsetfunc)0,                /* tp_descr_set */
        0,                              /* tp_dictoffset */
        (initproc)PP_ExClass1_tp_init,                   /* tp_init */
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

}  // namespace example
}  // namespace nested
