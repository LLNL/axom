// pyExClass2type.cpp
// This is generated code, do not edit
// blah blah
// yada yada
//
#include "pyUserLibrarymodule.hpp"
// splicer begin class.ExClass2.impl.include
// splicer end class.ExClass2.impl.include

namespace example {
namespace nested {
// splicer begin class.ExClass2.impl.C_definition
// splicer end class.ExClass2.impl.C_definition
// splicer begin class.ExClass2.impl.additional_methods
// splicer end class.ExClass2.impl.additional_methods
static void
PP_ExClass2_tp_dealloc (PP_ExClass2 *self)
{
// splicer begin class.ExClass2.type.dealloc
    PyErr_SetString(PyExc_NotImplementedError, "dealloc");
    return ;
// splicer end class.ExClass2.type.dealloc
}
static int
PP_ExClass2_tp_print (PP_ExClass2 *self, FILE *fp, int flags)
{
// splicer begin class.ExClass2.type.print
    PyErr_SetString(PyExc_NotImplementedError, "print");
    return -1;
// splicer end class.ExClass2.type.print
}
static int
PP_ExClass2_tp_compare (PP_ExClass2 *self, PyObject *)
{
// splicer begin class.ExClass2.type.compare
    PyErr_SetString(PyExc_NotImplementedError, "compare");
    return -1;
// splicer end class.ExClass2.type.compare
}
static PyObject *
PP_ExClass2_tp_getattr (PP_ExClass2 *self, char *name)
{
// splicer begin class.ExClass2.type.getattr
    PyErr_SetString(PyExc_NotImplementedError, "getattr");
    return NULL;
// splicer end class.ExClass2.type.getattr
}
static int
PP_ExClass2_tp_setattr (PP_ExClass2 *self, char *name, PyObject *value)
{
// splicer begin class.ExClass2.type.setattr
    PyErr_SetString(PyExc_NotImplementedError, "setattr");
    return -1;
// splicer end class.ExClass2.type.setattr
}
static PyObject *
PP_ExClass2_tp_getattro (PP_ExClass2 *self, PyObject *name)
{
// splicer begin class.ExClass2.type.getattro
    PyErr_SetString(PyExc_NotImplementedError, "getattro");
    return NULL;
// splicer end class.ExClass2.type.getattro
}
static int
PP_ExClass2_tp_setattro (PP_ExClass2 *self, PyObject *name, PyObject *value)
{
// splicer begin class.ExClass2.type.setattro
    PyErr_SetString(PyExc_NotImplementedError, "setattro");
    return -1;
// splicer end class.ExClass2.type.setattro
}
static PyObject *
PP_ExClass2_tp_repr (PP_ExClass2 *self)
{
// splicer begin class.ExClass2.type.repr
    PyErr_SetString(PyExc_NotImplementedError, "repr");
    return NULL;
// splicer end class.ExClass2.type.repr
}
static long
PP_ExClass2_tp_hash (PP_ExClass2 *self)
{
// splicer begin class.ExClass2.type.hash
    PyErr_SetString(PyExc_NotImplementedError, "hash");
    return -1;
// splicer end class.ExClass2.type.hash
}
static PyObject *
PP_ExClass2_tp_call (PP_ExClass2 *self, PyObject *args, PyObject *kwds)
{
// splicer begin class.ExClass2.type.call
    PyErr_SetString(PyExc_NotImplementedError, "call");
    return NULL;
// splicer end class.ExClass2.type.call
}
static PyObject *
PP_ExClass2_tp_str (PP_ExClass2 *self)
{
// splicer begin class.ExClass2.type.str
    PyErr_SetString(PyExc_NotImplementedError, "str");
    return NULL;
// splicer end class.ExClass2.type.str
}
static int
PP_ExClass2_tp_init (PP_ExClass2 *self, PyObject *args, PyObject *kwds)
{
// splicer begin class.ExClass2.type.init
    PyErr_SetString(PyExc_NotImplementedError, "init");
    return -1;
// splicer end class.ExClass2.type.init
}
static PyObject *
PP_ExClass2_tp_alloc (PyTypeObject *type, Py_ssize_t nitems)
{
// splicer begin class.ExClass2.type.alloc
    PyErr_SetString(PyExc_NotImplementedError, "alloc");
    return NULL;
// splicer end class.ExClass2.type.alloc
}
static PyObject *
PP_ExClass2_tp_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
// splicer begin class.ExClass2.type.new
    PyErr_SetString(PyExc_NotImplementedError, "new");
    return NULL;
// splicer end class.ExClass2.type.new
}
static void
PP_ExClass2_tp_free (void *op)
{
// splicer begin class.ExClass2.type.free
    PyErr_SetString(PyExc_NotImplementedError, "free");
    return ;
// splicer end class.ExClass2.type.free
}
static void
PP_ExClass2_tp_del (PP_ExClass2 *self)
{
// splicer begin class.ExClass2.type.del
    PyErr_SetString(PyExc_NotImplementedError, "del");
    return ;
// splicer end class.ExClass2.type.del
}

static char PP_exclass2_delete__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_delete(
  PP_ExClass2 *self,
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin class.ExClass2.method.delete
    delete self->BBB;
    self->BBB = NULL;
    Py_RETURN_NONE;
// splicer end class.ExClass2.method.delete
}

static char PP_exclass2_get_name__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_get_name(
  PP_ExClass2 *self,
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin class.ExClass2.method.get_name
    const std::string & rv = self->BBB->getName();
    PyObject * SH_Py_rv = PyString_FromString(rv.c_str());
    return (PyObject *) SH_Py_rv;
// splicer end class.ExClass2.method.get_name
}

static char PP_exclass2_get_name2__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_get_name2(
  PP_ExClass2 *self,
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin class.ExClass2.method.get_name2
    const std::string & rv = self->BBB->getName2();
    PyObject * SH_Py_rv = PyString_FromString(rv.c_str());
    return (PyObject *) SH_Py_rv;
// splicer end class.ExClass2.method.get_name2
}

static char PP_exclass2_get_name3__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_get_name3(
  PP_ExClass2 *self,
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin class.ExClass2.method.get_name3
    std::string & rv = self->BBB->getName3();
    PyObject * SH_Py_rv = PyString_FromString(rv.c_str());
    return (PyObject *) SH_Py_rv;
// splicer end class.ExClass2.method.get_name3
}

static char PP_exclass2_get_name4__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_get_name4(
  PP_ExClass2 *self,
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin class.ExClass2.method.get_name4
    std::string & rv = self->BBB->getName4();
    PyObject * SH_Py_rv = PyString_FromString(rv.c_str());
    return (PyObject *) SH_Py_rv;
// splicer end class.ExClass2.method.get_name4
}

static char PP_exclass2_get_name_length__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_get_name_length(
  PP_ExClass2 *self,
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin class.ExClass2.method.get_name_length
    const int rv = self->BBB->GetNameLength();
    return Py_BuildValue("i", rv);
// splicer end class.ExClass2.method.get_name_length
}

static char PP_exclass2_get_class1__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_get_class1(
  PP_ExClass2 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass2.method.get_class1
    const ExClass1 * in;
    PP_ExClass1 * SH_Py_in;
    const char *SH_kwcpp = "in";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!:get_class1", SH_kw_list,
        &PP_ExClass1_Type, &SH_Py_in))
    {
        return NULL;
    }
    in = SH_Py_in ? SH_Py_in->BBB : NULL;
    ExClass1 * rv = self->BBB->get_class1(in);
    PP_ExClass1 * SH_Py_rv = PyObject_New(PP_ExClass1, &PP_ExClass1_Type);
    SH_Py_rv->BBB = rv;
    return (PyObject *) SH_Py_rv;
// splicer end class.ExClass2.method.get_class1
}

static char PP_exclass2_declare_1__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_declare_1(
  PP_ExClass2 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass2.method.declare
    Py_ssize_t SH_nargs = 0;
    TypeID type;
    SidreLength len;
    const char *SH_kwcpp = "type\0len";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0,(char *) SH_kwcpp+5, NULL };

    if (args != NULL) SH_nargs += PyTuple_Size(args);
    if (kwds != NULL) SH_nargs += PyDict_Size(args);
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|l:declare", SH_kw_list,
        &type, &len))
    {
        return NULL;
    }
    switch (SH_nargs) {
    case 1:
        self->BBB->declare(getTypeID(type));
        break;
    case 2:
        self->BBB->declare(getTypeID(type), len);
        break;
    }
    Py_RETURN_NONE;
// splicer end class.ExClass2.method.declare
}

static char PP_exclass2_destroyall__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_destroyall(
  PP_ExClass2 *self,
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin class.ExClass2.method.destroyall
    self->BBB->destroyall();
    Py_RETURN_NONE;
// splicer end class.ExClass2.method.destroyall
}

static char PP_exclass2_get_type_id__doc__[] =
"documentation"
;

static PyObject *
PP_exclass2_get_type_id(
  PP_ExClass2 *self,
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin class.ExClass2.method.get_type_id
    TypeID rv = self->BBB->getTypeID();
    return Py_BuildValue("i", static_cast<int>(rv));
// splicer end class.ExClass2.method.get_type_id
}
// splicer begin class.ExClass2.impl.after_methods
// splicer end class.ExClass2.impl.after_methods
static PyMethodDef PP_ExClass2_methods[] = {
{"delete", (PyCFunction)PP_exclass2_delete, METH_NOARGS, PP_exclass2_delete__doc__},
{"getName", (PyCFunction)PP_exclass2_get_name, METH_NOARGS, PP_exclass2_get_name__doc__},
{"getName2", (PyCFunction)PP_exclass2_get_name2, METH_NOARGS, PP_exclass2_get_name2__doc__},
{"getName3", (PyCFunction)PP_exclass2_get_name3, METH_NOARGS, PP_exclass2_get_name3__doc__},
{"getName4", (PyCFunction)PP_exclass2_get_name4, METH_NOARGS, PP_exclass2_get_name4__doc__},
{"GetNameLength", (PyCFunction)PP_exclass2_get_name_length, METH_NOARGS, PP_exclass2_get_name_length__doc__},
{"get_class1", (PyCFunction)PP_exclass2_get_class1, METH_VARARGS|METH_KEYWORDS, PP_exclass2_get_class1__doc__},
{"declare", (PyCFunction)PP_exclass2_declare_1, METH_VARARGS|METH_KEYWORDS, PP_exclass2_declare_1__doc__},
{"destroyall", (PyCFunction)PP_exclass2_destroyall, METH_NOARGS, PP_exclass2_destroyall__doc__},
{"getTypeID", (PyCFunction)PP_exclass2_get_type_id, METH_NOARGS, PP_exclass2_get_type_id__doc__},
// splicer begin class.ExClass2.PyMethodDef
// splicer end class.ExClass2.PyMethodDef
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

static char ExClass2__doc__[] =
"virtual class"
;

/* static */
PyTypeObject PP_ExClass2_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "userlibrary.ExClass2",                       /* tp_name */
        sizeof(PP_ExClass2),         /* tp_basicsize */
        0,                              /* tp_itemsize */
        /* Methods to implement standard operations */
        (destructor)PP_ExClass2_tp_dealloc,                 /* tp_dealloc */
        (printfunc)PP_ExClass2_tp_print,                   /* tp_print */
        (getattrfunc)PP_ExClass2_tp_getattr,                 /* tp_getattr */
        (setattrfunc)PP_ExClass2_tp_setattr,                 /* tp_setattr */
#ifdef IS_PY3K
        0,                               /* tp_reserved */
#else
        (cmpfunc)PP_ExClass2_tp_compare,                     /* tp_compare */
#endif
        (reprfunc)PP_ExClass2_tp_repr,                    /* tp_repr */
        /* Method suites for standard classes */
        0,                              /* tp_as_number */
        0,                              /* tp_as_sequence */
        0,                              /* tp_as_mapping */
        /* More standard operations (here for binary compatibility) */
        (hashfunc)PP_ExClass2_tp_hash,                    /* tp_hash */
        (ternaryfunc)PP_ExClass2_tp_call,                 /* tp_call */
        (reprfunc)PP_ExClass2_tp_str,                    /* tp_str */
        (getattrofunc)PP_ExClass2_tp_getattro,                /* tp_getattro */
        (setattrofunc)PP_ExClass2_tp_setattro,                /* tp_setattro */
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
        PP_ExClass2_methods,                             /* tp_methods */
        0,                              /* tp_members */
        0,                             /* tp_getset */
        0,                              /* tp_base */
        0,                              /* tp_dict */
        (descrgetfunc)0,                /* tp_descr_get */
        (descrsetfunc)0,                /* tp_descr_set */
        0,                              /* tp_dictoffset */
        (initproc)PP_ExClass2_tp_init,                   /* tp_init */
        (allocfunc)PP_ExClass2_tp_alloc,                  /* tp_alloc */
        (newfunc)PP_ExClass2_tp_new,                    /* tp_new */
        (freefunc)PP_ExClass2_tp_free,                   /* tp_free */
        (inquiry)0,                     /* tp_is_gc */
        0,                              /* tp_bases */
        0,                              /* tp_mro */
        0,                              /* tp_cache */
        0,                              /* tp_subclasses */
        0,                              /* tp_weaklist */
        (destructor)PP_ExClass2_tp_del,                 /* tp_del */
        0,                              /* tp_version_tag */
#ifdef IS_PY3K
        (destructor)0,                  /* tp_finalize */
#endif
};

}  // namespace nested
}  // namespace example
