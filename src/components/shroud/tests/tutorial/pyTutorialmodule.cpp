// pyTutorialmodule.cpp
// This is generated code, do not edit
#include "pyTutorialmodule.hpp"
// splicer begin include
// splicer end include

namespace tutorial {
// splicer begin C_definition
// splicer end C_definition
PyObject *PY_error_obj;
// splicer begin additional_functions
// splicer end additional_functions

static char PY_function1__doc__[] =
"documentation"
;

static PyObject *
PY_function1(
  PyObject *,  // self unused
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin function.function1
    Function1();
    Py_RETURN_NONE;
// splicer end function.function1
}

static char PY_function2__doc__[] =
"documentation"
;

static PyObject *
PY_function2(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function2
    double arg1;
    int arg2;
    const char *SH_kwcpp = "arg1\0arg2";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0,(char *) SH_kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "di:Function2", SH_kw_list,
        &arg1, &arg2))
    {
        return NULL;
    }
    double rv = Function2(arg1, arg2);
    return Py_BuildValue("d", rv);
// splicer end function.function2
}

static char PY_function3__doc__[] =
"documentation"
;

static PyObject *
PY_function3(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function3
    bool arg;
    PyObject * SH_Py_arg;
    const char *SH_kwcpp = "arg";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!:Function3", SH_kw_list,
        &PyBool_Type, &SH_Py_arg))
    {
        return NULL;
    }
    arg = PyObject_IsTrue(SH_Py_arg);
    bool rv = Function3(arg);
    PyObject * SH_Py_rv = PyBool_FromLong(rv);
    return (PyObject *) SH_Py_rv;
// splicer end function.function3
}

static char PY_function4a__doc__[] =
"documentation"
;

static PyObject *
PY_function4a(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function4a
    const char * arg1;
    const char * arg2;
    const char *SH_kwcpp = "arg1\0arg2";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0,(char *) SH_kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ss:Function4a", SH_kw_list,
        &arg1, &arg2))
    {
        return NULL;
    }
    const std::string SH_arg1(arg1);
    const std::string SH_arg2(arg2);
    const std::string rv = Function4a(SH_arg1, SH_arg2);
    PyObject * SH_Py_rv = PyString_FromString(rv.c_str());
    return (PyObject *) SH_Py_rv;
// splicer end function.function4a
}

static char PY_function4b__doc__[] =
"documentation"
;

static PyObject *
PY_function4b(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function4b
    const char * arg1;
    const char * arg2;
    const char *SH_kwcpp = "arg1\0arg2";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0,(char *) SH_kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ss:Function4b", SH_kw_list,
        &arg1, &arg2))
    {
        return NULL;
    }
    const std::string SH_arg1(arg1);
    const std::string SH_arg2(arg2);
    const std::string & rv = Function4b(SH_arg1, SH_arg2);
    PyObject * SH_Py_rv = PyString_FromString(rv.c_str());
    return (PyObject *) SH_Py_rv;
// splicer end function.function4b
}

static char PY_function5_arg1_arg2__doc__[] =
"documentation"
;

static PyObject *
PY_function5_arg1_arg2(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function5
    Py_ssize_t SH_nargs = 0;
    double arg1;
    bool arg2;
    PyObject * SH_Py_arg2;
    const char *SH_kwcpp = "arg1\0arg2";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0,(char *) SH_kwcpp+5, NULL };
    double rv;
    
    if (args != NULL) SH_nargs += PyTuple_Size(args);
    if (kwds != NULL) SH_nargs += PyDict_Size(args);
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dO!:Function5", SH_kw_list,
        &arg1, &PyBool_Type, &SH_Py_arg2))
    {
        return NULL;
    }
    switch (SH_nargs) {
    case 0:
        rv = Function5();
        break;
    case 1:
        rv = Function5(arg1);
        break;
    case 2:
        arg2 = PyObject_IsTrue(SH_Py_arg2);
        rv = Function5(arg1, arg2);
        break;
    }
    return Py_BuildValue("d", rv);
// splicer end function.function5
}

static PyObject *
PY_function6_from_name(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function6_from_name
    const char * name;
    const char *SH_kwcpp = "name";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:Function6", SH_kw_list,
        &name))
    {
        return NULL;
    }
    const std::string SH_name(name);
    Function6(SH_name);
    Py_RETURN_NONE;
// splicer end function.function6_from_name
}

static PyObject *
PY_function6_from_index(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function6_from_index
    int indx;
    const char *SH_kwcpp = "indx";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:Function6", SH_kw_list,
        &indx))
    {
        return NULL;
    }
    Function6(indx);
    Py_RETURN_NONE;
// splicer end function.function6_from_index
}

static char PY_function9__doc__[] =
"documentation"
;

static PyObject *
PY_function9(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function9
    double arg;
    const char *SH_kwcpp = "arg";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "d:Function9", SH_kw_list,
        &arg))
    {
        return NULL;
    }
    Function9(arg);
    Py_RETURN_NONE;
// splicer end function.function9
}

static PyObject *
PY_function10_0(
  PyObject *,  // self unused
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin function.function10_0
    Function10();
    Py_RETURN_NONE;
// splicer end function.function10_0
}

static PyObject *
PY_function10_1(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function10_1
    const char * name;
    double arg2;
    const char *SH_kwcpp = "name\0arg2";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0,(char *) SH_kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "sd:Function10", SH_kw_list,
        &name, &arg2))
    {
        return NULL;
    }
    const std::string SH_name(name);
    Function10(SH_name, arg2);
    Py_RETURN_NONE;
// splicer end function.function10_1
}

static PyObject *
PY_overload1_num_offset_stride(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.overload1_num_offset_stride
    Py_ssize_t SH_nargs = 0;
    int num;
    int offset;
    int stride;
    const char *SH_kwcpp = "num\0offset\0stride";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0,(char *) SH_kwcpp+4,(char *) SH_kwcpp+11, NULL };
    int rv;
    
    if (args != NULL) SH_nargs += PyTuple_Size(args);
    if (kwds != NULL) SH_nargs += PyDict_Size(args);
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|ii:overload1", SH_kw_list,
        &num, &offset, &stride))
    {
        return NULL;
    }
    switch (SH_nargs) {
    case 1:
        rv = overload1(num);
        break;
    case 2:
        rv = overload1(num, offset);
        break;
    case 3:
        rv = overload1(num, offset, stride);
        break;
    }
    return Py_BuildValue("i", rv);
// splicer end function.overload1_num_offset_stride
}

static PyObject *
PY_overload1_5(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.overload1_5
    Py_ssize_t SH_nargs = 0;
    double type;
    int num;
    int offset;
    int stride;
    const char *SH_kwcpp = "type\0num\0offset\0stride";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0,(char *) SH_kwcpp+5,(char *) SH_kwcpp+9,(char *) SH_kwcpp+16, NULL };
    int rv;
    
    if (args != NULL) SH_nargs += PyTuple_Size(args);
    if (kwds != NULL) SH_nargs += PyDict_Size(args);
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "di|ii:overload1", SH_kw_list,
        &type, &num, &offset, &stride))
    {
        return NULL;
    }
    switch (SH_nargs) {
    case 2:
        rv = overload1(type, num);
        break;
    case 3:
        rv = overload1(type, num, offset);
        break;
    case 4:
        rv = overload1(type, num, offset, stride);
        break;
    }
    return Py_BuildValue("i", rv);
// splicer end function.overload1_5
}

static char PY_typefunc__doc__[] =
"documentation"
;

static PyObject *
PY_typefunc(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.typefunc
    TypeID arg;
    const char *SH_kwcpp = "arg";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:typefunc", SH_kw_list,
        &arg))
    {
        return NULL;
    }
    TypeID rv = typefunc(arg);
    return Py_BuildValue("i", rv);
// splicer end function.typefunc
}

static char PY_enumfunc__doc__[] =
"documentation"
;

static PyObject *
PY_enumfunc(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.enumfunc
    EnumTypeID arg;
    const char *SH_kwcpp = "arg";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:enumfunc", SH_kw_list,
        &arg))
    {
        return NULL;
    }
    EnumTypeID rv = enumfunc(static_cast<EnumTypeID>(arg));
    return Py_BuildValue("i", static_cast<int>(rv));
// splicer end function.enumfunc
}

static char PY_useclass__doc__[] =
"documentation"
;

static PyObject *
PY_useclass(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.useclass
    const Class1 * arg1;
    PY_Class1 * SH_Py_arg1;
    const char *SH_kwcpp = "arg1";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!:useclass", SH_kw_list,
        &PY_Class1_Type, &SH_Py_arg1))
    {
        return NULL;
    }
    arg1 = SH_Py_arg1 ? SH_Py_arg1->BBB : NULL;
    useclass(arg1);
    Py_RETURN_NONE;
// splicer end function.useclass
}

static char PY_last_function_called__doc__[] =
"documentation"
;

static PyObject *
PY_last_function_called(
  PyObject *,  // self unused
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin function.last_function_called
    const std::string & rv = LastFunctionCalled();
    PyObject * SH_Py_rv = PyString_FromString(rv.c_str());
    return (PyObject *) SH_Py_rv;
// splicer end function.last_function_called
}

static char PY_function10__doc__[] =
"documentation"
;

static PyObject *
PY_function10(
  PyObject *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function10
    Py_ssize_t SH_nargs = 0;
    if (args != NULL) SH_nargs += PyTuple_Size(args);
    if (kwds != NULL) SH_nargs += PyDict_Size(args);
    PyObject *rvobj;
    if (SH_nargs == 0) {
        rvobj = PY_function10_0(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    if (SH_nargs == 2) {
        rvobj = PY_function10_1(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
    return NULL;
// splicer end function.function10
}

static char PY_function6__doc__[] =
"documentation"
;

static PyObject *
PY_function6(
  PyObject *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function6
    Py_ssize_t SH_nargs = 0;
    if (args != NULL) SH_nargs += PyTuple_Size(args);
    if (kwds != NULL) SH_nargs += PyDict_Size(args);
    PyObject *rvobj;
    if (SH_nargs == 1) {
        rvobj = PY_function6_from_name(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    if (SH_nargs == 1) {
        rvobj = PY_function6_from_index(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
    return NULL;
// splicer end function.function6
}

static char PY_overload1__doc__[] =
"documentation"
;

static PyObject *
PY_overload1(
  PyObject *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.overload1
    Py_ssize_t SH_nargs = 0;
    if (args != NULL) SH_nargs += PyTuple_Size(args);
    if (kwds != NULL) SH_nargs += PyDict_Size(args);
    PyObject *rvobj;
    if (SH_nargs >= 1 && SH_nargs <= 3) {
        rvobj = PY_overload1_num_offset_stride(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    if (SH_nargs >= 2 && SH_nargs <= 4) {
        rvobj = PY_overload1_5(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
    return NULL;
// splicer end function.overload1
}
static PyMethodDef PY_methods[] = {
{"Function1", (PyCFunction)PY_function1, METH_NOARGS, PY_function1__doc__},
{"Function2", (PyCFunction)PY_function2, METH_VARARGS|METH_KEYWORDS, PY_function2__doc__},
{"Function3", (PyCFunction)PY_function3, METH_VARARGS|METH_KEYWORDS, PY_function3__doc__},
{"Function4a", (PyCFunction)PY_function4a, METH_VARARGS|METH_KEYWORDS, PY_function4a__doc__},
{"Function4b", (PyCFunction)PY_function4b, METH_VARARGS|METH_KEYWORDS, PY_function4b__doc__},
{"Function5", (PyCFunction)PY_function5_arg1_arg2, METH_VARARGS|METH_KEYWORDS, PY_function5_arg1_arg2__doc__},
{"Function9", (PyCFunction)PY_function9, METH_VARARGS|METH_KEYWORDS, PY_function9__doc__},
{"typefunc", (PyCFunction)PY_typefunc, METH_VARARGS|METH_KEYWORDS, PY_typefunc__doc__},
{"enumfunc", (PyCFunction)PY_enumfunc, METH_VARARGS|METH_KEYWORDS, PY_enumfunc__doc__},
{"useclass", (PyCFunction)PY_useclass, METH_VARARGS|METH_KEYWORDS, PY_useclass__doc__},
{"LastFunctionCalled", (PyCFunction)PY_last_function_called, METH_NOARGS, PY_last_function_called__doc__},
{"Function10", (PyCFunction)PY_function10, METH_VARARGS|METH_KEYWORDS, PY_function10__doc__},
{"Function6", (PyCFunction)PY_function6, METH_VARARGS|METH_KEYWORDS, PY_function6__doc__},
{"overload1", (PyCFunction)PY_overload1, METH_VARARGS|METH_KEYWORDS, PY_overload1__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

/*
 * inittutorial - Initialization function for the module
 * *must* be called inittutorial
 */
static char PY__doc__[] =
"library documentation"
;

struct module_state {
    PyObject *error;
};

#ifdef IS_PY3K
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

#ifdef IS_PY3K
static int tutorial_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int tutorial_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "tutorial", /* m_name */
    PY__doc__, /* m_doc */
    sizeof(struct module_state), /* m_size */
    PY_methods, /* m_methods */
    NULL, /* m_reload */
    tutorial_traverse, /* m_traverse */
    tutorial_clear, /* m_clear */
    NULL  /* m_free */
};

#define RETVAL m
#define INITERROR return NULL
#else
#define RETVAL
#define INITERROR return
#endif

#ifdef __cplusplus
extern "C" {
#endif
PyMODINIT_FUNC
MOD_INITBASIS(void)
{
    PyObject *m = NULL;
    const char * error_name = "tutorial.Error";

// splicer begin C_init_locals
// splicer end C_init_locals


    /* Create the module and add the functions */
#ifdef IS_PY3K
    m = PyModule_Create(&moduledef);
#else
    m = Py_InitModule4("tutorial", PY_methods,
                       PY__doc__,
                       (PyObject*)NULL,PYTHON_API_VERSION);
#endif
    if (m == NULL)
        return RETVAL;
    struct module_state *st = GETSTATE(m);


// Class1
    PY_Class1_Type.tp_new   = PyType_GenericNew;
    PY_Class1_Type.tp_alloc = PyType_GenericAlloc;
    if (PyType_Ready(&PY_Class1_Type) < 0)
        return RETVAL;
    Py_INCREF(&PY_Class1_Type);
    PyModule_AddObject(m, "Class1", (PyObject *)&PY_Class1_Type);


    PY_error_obj = PyErr_NewException((char *) error_name, NULL, NULL);
    if (PY_error_obj == NULL)
        return RETVAL;
    st->error = PY_error_obj;
    PyModule_AddObject(m, "Error", st->error);

// splicer begin C_init_body
// splicer end C_init_body

    /* Check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module tutorial");
    return RETVAL;
}
#ifdef __cplusplus
}
#endif

}  // namespace tutorial
