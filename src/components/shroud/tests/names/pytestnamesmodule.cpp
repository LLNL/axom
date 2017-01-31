// pytestnamesmodule.cpp
// This is generated code, do not edit
#include "pytestnamesmodule.hpp"
// splicer begin include
// splicer end include

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
    function1();
    Py_RETURN_NONE;
// splicer end function.function1
}

static char PY_function2__doc__[] =
"documentation"
;

static PyObject *
PY_function2(
  PyObject *,  // self unused
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin function.function2
    function2();
    Py_RETURN_NONE;
// splicer end function.function2
}

static PyObject *
PY_function3a_0(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function3a_0
    int i;
    const char *SH_kwcpp = "i";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:function3a", SH_kw_list,
        &i))
    {
        return NULL;
    }
    function3a(i);
    Py_RETURN_NONE;
// splicer end function.function3a_0
}

static PyObject *
PY_function3a_1(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function3a_1
    long i;
    const char *SH_kwcpp = "i";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "l:function3a", SH_kw_list,
        &i))
    {
        return NULL;
    }
    function3a(i);
    Py_RETURN_NONE;
// splicer end function.function3a_1
}

static char PY_function3a__doc__[] =
"documentation"
;

static PyObject *
PY_function3a(
  PyObject *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.function3a
    Py_ssize_t SH_nargs = 0;
    if (args != NULL) SH_nargs += PyTuple_Size(args);
    if (kwds != NULL) SH_nargs += PyDict_Size(args);
    PyObject *rvobj;
    if (SH_nargs == 1) {
        rvobj = PY_function3a_0(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    if (SH_nargs == 1) {
        rvobj = PY_function3a_1(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
    return NULL;
// splicer end function.function3a
}
static PyMethodDef PY_methods[] = {
{"function1", (PyCFunction)PY_function1, METH_NOARGS, PY_function1__doc__},
{"function2", (PyCFunction)PY_function2, METH_NOARGS, PY_function2__doc__},
{"function3a", (PyCFunction)PY_function3a, METH_VARARGS|METH_KEYWORDS, PY_function3a__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

/*
 * inittestnames - Initialization function for the module
 * *must* be called inittestnames
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
static int testnames_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int testnames_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "testnames", /* m_name */
    PY__doc__, /* m_doc */
    sizeof(struct module_state), /* m_size */
    PY_methods, /* m_methods */
    NULL, /* m_reload */
    testnames_traverse, /* m_traverse */
    testnames_clear, /* m_clear */
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
    const char * error_name = "testnames.Error";

// splicer begin C_init_locals
// splicer end C_init_locals


    /* Create the module and add the functions */
#ifdef IS_PY3K
    m = PyModule_Create(&moduledef);
#else
    m = Py_InitModule4("testnames", PY_methods,
                       PY__doc__,
                       (PyObject*)NULL,PYTHON_API_VERSION);
#endif
    if (m == NULL)
        return RETVAL;
    struct module_state *st = GETSTATE(m);


// Names
    PY_Names_Type.tp_new   = PyType_GenericNew;
    PY_Names_Type.tp_alloc = PyType_GenericAlloc;
    if (PyType_Ready(&PY_Names_Type) < 0)
        return RETVAL;
    Py_INCREF(&PY_Names_Type);
    PyModule_AddObject(m, "Names", (PyObject *)&PY_Names_Type);


    PY_error_obj = PyErr_NewException((char *) error_name, NULL, NULL);
    if (PY_error_obj == NULL)
        return RETVAL;
    st->error = PY_error_obj;
    PyModule_AddObject(m, "Error", st->error);

// splicer begin C_init_body
// splicer end C_init_body

    /* Check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module testnames");
    return RETVAL;
}
#ifdef __cplusplus
}
#endif

