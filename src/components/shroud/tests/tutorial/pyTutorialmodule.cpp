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
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
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
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function2
    double arg1;
    int arg2;
    const char *kwcpp = "arg1\0arg2";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "di:Function2", kw_list,
        &arg1, &arg2))
    {
        return NULL;
    }
    double rv = Function2(arg1, arg2);
    return Py_BuildValue("d", rv);
// splicer end function.function2
}

static char PY_sum__doc__[] =
"documentation"
;

static PyObject *
PY_sum(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.sum
    int len;
    int * values;
    int * result;
    const char *kwcpp = "len\0values\0result";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+4,(char *) kwcpp+11, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "iii:Sum", kw_list,
        &len, &values, &result))
    {
        return NULL;
    }
    Sum(len, values, result);
    Py_RETURN_NONE;
// splicer end function.sum
}

static char PY_function3__doc__[] =
"documentation"
;

static PyObject *
PY_function3(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function3
    bool arg;
    const char *kwcpp = "arg";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O:Function3", kw_list,
        &arg))
    {
        return NULL;
    }
    bool rv = Function3(arg);
    return PyBool_FromLong(rv);
// splicer end function.function3
}

static char PY_function4a__doc__[] =
"documentation"
;

static PyObject *
PY_function4a(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function4a
    const char * arg1;
    const char * arg2;
    const char *kwcpp = "arg1\0arg2";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ss:Function4a", kw_list,
        &arg1, &arg2))
    {
        return NULL;
    }
    const std::string & rv = Function4a(arg1, arg2);
    return PyString_FromString(rv.c_str());
// splicer end function.function4a
}

static char PY_function4b__doc__[] =
"documentation"
;

static PyObject *
PY_function4b(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function4b
    const char * arg1;
    const char * arg2;
    const char *kwcpp = "arg1\0arg2";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ss:Function4b", kw_list,
        &arg1, &arg2))
    {
        return NULL;
    }
    const std::string & rv = Function4b(arg1, arg2);
    return PyString_FromString(rv.c_str());
// splicer end function.function4b
}

static char PY_function5_arg1_arg2__doc__[] =
"documentation"
;

static PyObject *
PY_function5_arg1_arg2(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function5_arg1_arg2
    double arg1;
    bool arg2;
    const char *kwcpp = "arg1\0arg2";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    arg1 = 3.1415;
    arg2 = true;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|dO:Function5", kw_list,
        &arg1, &arg2))
    {
        return NULL;
    }
    double rv = Function5(arg1, arg2);
    return Py_BuildValue("d", rv);
// splicer end function.function5_arg1_arg2
}

static char PY_function6_from_name__doc__[] =
"documentation"
;

static PyObject *
PY_function6_from_name(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function6_from_name
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:Function6", kw_list,
        &name))
    {
        return NULL;
    }
    Function6(name);
    Py_RETURN_NONE;
// splicer end function.function6_from_name
}

static char PY_function6_from_index__doc__[] =
"documentation"
;

static PyObject *
PY_function6_from_index(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function6_from_index
    int indx;
    const char *kwcpp = "indx";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:Function6", kw_list,
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
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function9
    double arg;
    const char *kwcpp = "arg";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "d:Function9", kw_list,
        &arg))
    {
        return NULL;
    }
    Function9(arg);
    Py_RETURN_NONE;
// splicer end function.function9
}

static char PY_function10_0__doc__[] =
"documentation"
;

static PyObject *
PY_function10_0(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function10_0
    Function10();
    Py_RETURN_NONE;
// splicer end function.function10_0
}

static char PY_function10_1__doc__[] =
"documentation"
;

static PyObject *
PY_function10_1(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function10_1
    const char * name;
    double arg2;
    const char *kwcpp = "name\0arg2";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "sd:Function10", kw_list,
        &name, &arg2))
    {
        return NULL;
    }
    Function10(name, arg2);
    Py_RETURN_NONE;
// splicer end function.function10_1
}

static char PY_overload1_num_offset_stride__doc__[] =
"documentation"
;

static PyObject *
PY_overload1_num_offset_stride(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.overload1_num_offset_stride
    int num;
    int offset;
    int stride;
    const char *kwcpp = "num\0offset\0stride";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+4,(char *) kwcpp+11, NULL };
    
    offset = 0;
    stride = 1;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i|ii:overload1", kw_list,
        &num, &offset, &stride))
    {
        return NULL;
    }
    int rv = overload1(num, offset, stride);
    return Py_BuildValue("i", rv);
// splicer end function.overload1_num_offset_stride
}

static char PY_overload1_5__doc__[] =
"documentation"
;

static PyObject *
PY_overload1_5(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.overload1_5
    double type;
    int num;
    int offset;
    int stride;
    const char *kwcpp = "type\0num\0offset\0stride";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5,(char *) kwcpp+9,(char *) kwcpp+16, NULL };
    
    offset = 0;
    stride = 1;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "di|ii:overload1", kw_list,
        &type, &num, &offset, &stride))
    {
        return NULL;
    }
    int rv = overload1(type, num, offset, stride);
    return Py_BuildValue("i", rv);
// splicer end function.overload1_5
}

static char PY_typefunc__doc__[] =
"documentation"
;

static PyObject *
PY_typefunc(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.typefunc
    int arg;
    const char *kwcpp = "arg";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:typefunc", kw_list,
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
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.enumfunc
    int arg;
    const char *kwcpp = "arg";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:enumfunc", kw_list,
        &arg))
    {
        return NULL;
    }
    EnumTypeID rv = enumfunc(static_cast<EnumTypeID>(arg));
    return Py_BuildValue("i", rv);
// splicer end function.enumfunc
}

static char PY_last_function_called__doc__[] =
"documentation"
;

static PyObject *
PY_last_function_called(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.last_function_called
    const std::string & rv = LastFunctionCalled();
    return PyString_FromString(rv.c_str());
// splicer end function.last_function_called
}
static PyMethodDef PY_methods[] = {
{"Function1", (PyCFunction)PY_function1, METH_NOARGS, PY_function1__doc__},
{"Function2", (PyCFunction)PY_function2, METH_VARARGS|METH_KEYWORDS, PY_function2__doc__},
{"Sum", (PyCFunction)PY_sum, METH_VARARGS|METH_KEYWORDS, PY_sum__doc__},
{"Function3", (PyCFunction)PY_function3, METH_VARARGS|METH_KEYWORDS, PY_function3__doc__},
{"Function4a", (PyCFunction)PY_function4a, METH_VARARGS|METH_KEYWORDS, PY_function4a__doc__},
{"Function4b", (PyCFunction)PY_function4b, METH_VARARGS|METH_KEYWORDS, PY_function4b__doc__},
{"Function5_arg1_arg2", (PyCFunction)PY_function5_arg1_arg2, METH_VARARGS|METH_KEYWORDS, PY_function5_arg1_arg2__doc__},
{"Function6_from_name", (PyCFunction)PY_function6_from_name, METH_VARARGS|METH_KEYWORDS, PY_function6_from_name__doc__},
{"Function6_from_index", (PyCFunction)PY_function6_from_index, METH_VARARGS|METH_KEYWORDS, PY_function6_from_index__doc__},
{"Function9", (PyCFunction)PY_function9, METH_VARARGS|METH_KEYWORDS, PY_function9__doc__},
{"Function10_0", (PyCFunction)PY_function10_0, METH_NOARGS, PY_function10_0__doc__},
{"Function10_1", (PyCFunction)PY_function10_1, METH_VARARGS|METH_KEYWORDS, PY_function10_1__doc__},
{"overload1_num_offset_stride", (PyCFunction)PY_overload1_num_offset_stride, METH_VARARGS|METH_KEYWORDS, PY_overload1_num_offset_stride__doc__},
{"overload1_5", (PyCFunction)PY_overload1_5, METH_VARARGS|METH_KEYWORDS, PY_overload1_5__doc__},
{"typefunc", (PyCFunction)PY_typefunc, METH_VARARGS|METH_KEYWORDS, PY_typefunc__doc__},
{"enumfunc", (PyCFunction)PY_enumfunc, METH_VARARGS|METH_KEYWORDS, PY_enumfunc__doc__},
{"LastFunctionCalled", (PyCFunction)PY_last_function_called, METH_NOARGS, PY_last_function_called__doc__},
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
