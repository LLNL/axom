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
    int incr;
    const char *kwcpp = "incr";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:incrementCount", kw_list,
        &incr))
    {
        return NULL;
    }
    int rv = self->BBB->incrementCount(incr);
    return Py_BuildValue("i", &rv);
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
    const std::string & rv = self->BBB->getName();
    if (! isNameValid(rv)) {
        PyErr_SetString(PyExc_KeyError, "'rv'");
        return NULL;
    }
    
    return Py_BuildValue("s", rv.c_str());
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
// splicer begin class.ExClass1.method.GetNameLength
    int rv = self->BBB->GetNameLength();
    return Py_BuildValue("i", &rv);
// splicer end class.ExClass1.method.GetNameLength
}

static char PP_exclass1_get_name_error_check__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_get_name_error_check(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.getNameErrorCheck
    const std::string & rv = self->BBB->getNameErrorCheck();
    return Py_BuildValue("s", rv.c_str());
// splicer end class.ExClass1.method.getNameErrorCheck
}

static char PP_exclass1_get_name_arg__doc__[] =
"documentation"
;

static PyObject *
PP_exclass1_get_name_arg(
  PP_ExClass1 *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.ExClass1.method.getNameArg
    const std::string & rv = self->BBB->getNameArg();
    return Py_BuildValue("s", rv.c_str());
// splicer end class.ExClass1.method.getNameArg
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
    ExClass2 * rv = self->BBB->getRoot();
    PP_ExClass2 * rv_obj = PyObject_New(PP_ExClass2, &PP_ExClass2_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
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
    int value;
    const char *kwcpp = "value";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:getValue", kw_list,
        &value))
    {
        return NULL;
    }
    int rv = self->BBB->getValue(value);
    return Py_BuildValue("i", &rv);
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
    long value;
    const char *kwcpp = "value";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "l:getValue", kw_list,
        &value))
    {
        return NULL;
    }
    long rv = self->BBB->getValue(value);
    return Py_BuildValue("l", &rv);
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
    void * rv = self->BBB->getAddr();
    return Py_BuildValue("O", rv);
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
    bool in;
    const char *kwcpp = "in";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O:hasAddr", kw_list,
        &in))
    {
        return NULL;
    }
    bool rv = self->BBB->hasAddr(in);
    return Py_BuildValue("O", &rv);
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
    self->BBB->SplicerSpecial();
    Py_RETURN_NONE;
// splicer end class.ExClass1.method.SplicerSpecial
}
static PyMethodDef PP_ExClass1_methods[] = {
{"incrementCount", (PyCFunction)PP_exclass1_increment_count, METH_VARARGS|METH_KEYWORDS, PP_exclass1_increment_count__doc__},
{"getName", (PyCFunction)PP_exclass1_get_name, METH_NOARGS, PP_exclass1_get_name__doc__},
{"GetNameLength", (PyCFunction)PP_exclass1_get_name_length, METH_NOARGS, PP_exclass1_get_name_length__doc__},
{"getNameErrorCheck", (PyCFunction)PP_exclass1_get_name_error_check, METH_NOARGS, PP_exclass1_get_name_error_check__doc__},
{"getNameArg", (PyCFunction)PP_exclass1_get_name_arg, METH_NOARGS, PP_exclass1_get_name_arg__doc__},
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
