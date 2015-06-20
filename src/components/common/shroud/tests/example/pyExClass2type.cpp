// pyExClass2type.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
#include "pyUserLibrarymodule.hpp"
// splicer begin class.C_definition
// splicer end class.C_definition
// splicer begin class.extra_methods
// splicer end class.extra_methods

static char PP_exclass2_ex_class2__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_ex_class2(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass2.method.ex_class2
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass2.method.ex_class2
}

static char PP_exclass2_ex_class1__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_ex_class1(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass2.method.ex_class1
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass2.method.ex_class1
}

static char PP_exclass2_get_name__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_get_name(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass2.method.get_name
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass2.method.get_name
}

static char PP_exclass2_get_name_length__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_get_name_length(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass2.method.get_name_length
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass2.method.get_name_length
}

static char PP_exclass2_get_class1__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_get_class1(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass2.method.get_class1
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass2.method.get_class1
}

static char PP_exclass2_declare__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_declare(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass2.method.declare
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass2.method.declare
}

static char PP_exclass2_destroyall__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_destroyall(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass2.method.destroyall
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass2.method.destroyall
}

static char PP_exclass2_get_type_id__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_get_type_id(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass2.method.get_type_id
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.ExClass2.method.get_type_id
}
static PyMethodDef PB_methods[] = {
{"ExClass2", (PyCFunction)PP_exclass2_ex_class2, METH_VARARGS|METH_KEYWORDS, PP_exclass2_ex_class2__doc__},
{"ExClass1", (PyCFunction)PP_exclass2_ex_class1, METH_VARARGS|METH_KEYWORDS, PP_exclass2_ex_class1__doc__},
{"getName", (PyCFunction)PP_exclass2_get_name, METH_VARARGS|METH_KEYWORDS, PP_exclass2_get_name__doc__},
{"get_name_length", (PyCFunction)PP_exclass2_get_name_length, METH_VARARGS|METH_KEYWORDS, PP_exclass2_get_name_length__doc__},
{"get_class1", (PyCFunction)PP_exclass2_get_class1, METH_VARARGS|METH_KEYWORDS, PP_exclass2_get_class1__doc__},
{"declare", (PyCFunction)PP_exclass2_declare, METH_VARARGS|METH_KEYWORDS, PP_exclass2_declare__doc__},
{"destroyall", (PyCFunction)PP_exclass2_destroyall, METH_VARARGS|METH_KEYWORDS, PP_exclass2_destroyall__doc__},
{"getTypeID", (PyCFunction)PP_exclass2_get_type_id, METH_VARARGS|METH_KEYWORDS, PP_exclass2_get_type_id__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

static char ExClass2__doc__[] =
"virtual class"
;

/* static */
PyTypeObject PP_ExClass2_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "basis.Dbnode",                       /* tp_name */
        sizeof(PP_ExClass2),         /* tp_basicsize */
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
        ExClass2__doc__,         /* tp_doc */
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
        0,                             /* tp_methods */
        0,                              /* tp_members */
        0,                             /* tp_getset */
        0,                              /* tp_base */
        0,                              /* tp_dict */
        (descrgetfunc)0,                /* tp_descr_get */
        (descrsetfunc)0,                /* tp_descr_set */
        0,                              /* tp_dictoffset */
        (initproc)0,                   /* tp_init */
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

