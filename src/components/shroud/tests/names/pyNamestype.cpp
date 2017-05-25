// pyNamestype.cpp
// This is generated code, do not edit
#include "pytestnamesmodule.hpp"
// splicer begin class.Names.impl.include
// splicer end class.Names.impl.include
// splicer begin class.Names.impl.C_definition
// splicer end class.Names.impl.C_definition
// splicer begin class.Names.impl.additional_methods
// splicer end class.Names.impl.additional_methods

static char PY_names_method1__doc__[] =
"documentation"
;

static PyObject *
PY_names_method1(
  PY_Names *self,
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin class.Names.method.method1
    self->BBB->method1();
    Py_RETURN_NONE;
// splicer end class.Names.method.method1
}

static char PY_names_method2__doc__[] =
"documentation"
;

static PyObject *
PY_names_method2(
  PY_Names *self,
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin class.Names.method.method2
    self->BBB->method2();
    Py_RETURN_NONE;
// splicer end class.Names.method.method2
}
// splicer begin class.Names.impl.after_methods
// splicer end class.Names.impl.after_methods
static PyMethodDef PY_Names_methods[] = {
{"method1", (PyCFunction)PY_names_method1, METH_NOARGS, PY_names_method1__doc__},
{"method2", (PyCFunction)PY_names_method2, METH_NOARGS, PY_names_method2__doc__},
// splicer begin class.Names.PyMethodDef
// splicer end class.Names.PyMethodDef
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

static char Names__doc__[] =
"virtual class"
;

/* static */
PyTypeObject PY_Names_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "testnames.Names",                       /* tp_name */
        sizeof(PY_Names),         /* tp_basicsize */
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
        Names__doc__,         /* tp_doc */
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
        PY_Names_methods,                             /* tp_methods */
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
