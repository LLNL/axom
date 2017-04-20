// pystringsmodule.cpp
// This is generated code, do not edit
#include "pystringsmodule.hpp"
// splicer begin include
// splicer end include

// splicer begin C_definition
// splicer end C_definition
PyObject *PY_error_obj;
// splicer begin additional_functions
// splicer end additional_functions

static char PY_pass_char__doc__[] =
"documentation"
;

static PyObject *
PY_pass_char(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.pass_char
    char status;
    const char *SH_kwcpp = "status";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:passChar", SH_kw_list,
        &status))
    {
        return NULL;
    }
    passChar(status);
    Py_RETURN_NONE;
// splicer end function.pass_char
}

static char PY_return_char__doc__[] =
"documentation"
;

static PyObject *
PY_return_char(
  PyObject *,  // self unused
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin function.return_char
    char rv = returnChar();
    PyObject * SH_Py_rv = PyString_FromString(rv);
    return (PyObject *) SH_Py_rv;
// splicer end function.return_char
}

static char PY_pass_char_ptr__doc__[] =
"documentation"
;

static PyObject *
PY_pass_char_ptr(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.pass_char_ptr
    char * dest;
    const char * src;
    const char *SH_kwcpp = "src";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:passCharPtr", SH_kw_list,
        &src))
    {
        return NULL;
    }
    passCharPtr(dest, src);
    PyObject * SH_Py_dest = PyString_FromString(dest);
    return (PyObject *) *SH_Py_dest;
// splicer end function.pass_char_ptr
}

static char PY_get_char1__doc__[] =
"documentation"
;

static PyObject *
PY_get_char1(
  PyObject *,  // self unused
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin function.get_char1
    const char * rv = getChar1();
    PyObject * SH_Py_rv = PyString_FromString(rv);
    return (PyObject *) SH_Py_rv;
// splicer end function.get_char1
}

static char PY_get_char2__doc__[] =
"documentation"
;

static PyObject *
PY_get_char2(
  PyObject *,  // self unused
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin function.get_char2
    const char * rv = getChar2();
    PyObject * SH_Py_rv = PyString_FromString(rv);
    return (PyObject *) SH_Py_rv;
// splicer end function.get_char2
}

static char PY_get_char3__doc__[] =
"documentation"
;

static PyObject *
PY_get_char3(
  PyObject *,  // self unused
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin function.get_char3
    const char * rv = getChar3();
    PyObject * SH_Py_rv = PyString_FromString(rv);
    return (PyObject *) SH_Py_rv;
// splicer end function.get_char3
}

static char PY_get_string1__doc__[] =
"documentation"
;

static PyObject *
PY_get_string1(
  PyObject *,  // self unused
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin function.get_string1
    const std::string & rv = getString1();
    PyObject * SH_Py_rv = PyString_FromString(rv.c_str());
    return (PyObject *) SH_Py_rv;
// splicer end function.get_string1
}

static char PY_get_string2__doc__[] =
"documentation"
;

static PyObject *
PY_get_string2(
  PyObject *,  // self unused
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin function.get_string2
    const std::string & rv = getString2();
    PyObject * SH_Py_rv = PyString_FromString(rv.c_str());
    return (PyObject *) SH_Py_rv;
// splicer end function.get_string2
}

static char PY_get_string3__doc__[] =
"documentation"
;

static PyObject *
PY_get_string3(
  PyObject *,  // self unused
  PyObject *,  // args unused
  PyObject *)  // kwds unused
{
// splicer begin function.get_string3
    const std::string & rv = getString3();
    PyObject * SH_Py_rv = PyString_FromString(rv.c_str());
    return (PyObject *) SH_Py_rv;
// splicer end function.get_string3
}

static char PY_accept_string_const_reference__doc__[] =
"documentation"
;

static PyObject *
PY_accept_string_const_reference(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.accept_string_const_reference
    const char * arg1;
    const char *SH_kwcpp = "arg1";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:acceptStringConstReference", SH_kw_list,
        &arg1))
    {
        return NULL;
    }
    const std::string SH_arg1(arg1);
    acceptStringConstReference(SH_arg1);
    Py_RETURN_NONE;
// splicer end function.accept_string_const_reference
}

static char PY_accept_string_reference__doc__[] =
"documentation"
;

static PyObject *
PY_accept_string_reference(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.accept_string_reference
    char * arg1;
    const char *SH_kwcpp = "arg1";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:acceptStringReference", SH_kw_list,
        &arg1))
    {
        return NULL;
    }
    std::string SH_arg1(arg1);
    acceptStringReference(SH_arg1);
    PyObject * SH_Py_arg1 = PyString_FromString(SH_arg1.c_str());
    return (PyObject *) *SH_Py_arg1;
// splicer end function.accept_string_reference
}

static char PY_explicit1__doc__[] =
"documentation"
;

static PyObject *
PY_explicit1(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.explicit1
    char * name;
    const char *SH_kwcpp = "name";
    char *SH_kw_list[] = { (char *) SH_kwcpp+0, NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:explicit1", SH_kw_list,
        &name))
    {
        return NULL;
    }
    explicit1(name);
    Py_RETURN_NONE;
// splicer end function.explicit1
}

static char PY_explicit2__doc__[] =
"documentation"
;

static PyObject *
PY_explicit2(
  PyObject *,  // self unused
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.explicit2
    char * name;
    const char *SH_kwcpp = "";
    char *SH_kw_list[] = { , NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, ":explicit2", SH_kw_list,
        ))
    {
        return NULL;
    }
    explicit2(name);
    PyObject * SH_Py_name = PyString_FromString(name);
    return (PyObject *) *SH_Py_name;
// splicer end function.explicit2
}
static PyMethodDef PY_methods[] = {
{"passChar", (PyCFunction)PY_pass_char, METH_VARARGS|METH_KEYWORDS, PY_pass_char__doc__},
{"returnChar", (PyCFunction)PY_return_char, METH_NOARGS, PY_return_char__doc__},
{"passCharPtr", (PyCFunction)PY_pass_char_ptr, METH_VARARGS|METH_KEYWORDS, PY_pass_char_ptr__doc__},
{"getChar1", (PyCFunction)PY_get_char1, METH_NOARGS, PY_get_char1__doc__},
{"getChar2", (PyCFunction)PY_get_char2, METH_NOARGS, PY_get_char2__doc__},
{"getChar3", (PyCFunction)PY_get_char3, METH_NOARGS, PY_get_char3__doc__},
{"getString1", (PyCFunction)PY_get_string1, METH_NOARGS, PY_get_string1__doc__},
{"getString2", (PyCFunction)PY_get_string2, METH_NOARGS, PY_get_string2__doc__},
{"getString3", (PyCFunction)PY_get_string3, METH_NOARGS, PY_get_string3__doc__},
{"acceptStringConstReference", (PyCFunction)PY_accept_string_const_reference, METH_VARARGS|METH_KEYWORDS, PY_accept_string_const_reference__doc__},
{"acceptStringReference", (PyCFunction)PY_accept_string_reference, METH_VARARGS|METH_KEYWORDS, PY_accept_string_reference__doc__},
{"explicit1", (PyCFunction)PY_explicit1, METH_VARARGS|METH_KEYWORDS, PY_explicit1__doc__},
{"explicit2", (PyCFunction)PY_explicit2, METH_VARARGS|METH_KEYWORDS, PY_explicit2__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

/*
 * initstrings - Initialization function for the module
 * *must* be called initstrings
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
static int strings_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int strings_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "strings", /* m_name */
    PY__doc__, /* m_doc */
    sizeof(struct module_state), /* m_size */
    PY_methods, /* m_methods */
    NULL, /* m_reload */
    strings_traverse, /* m_traverse */
    strings_clear, /* m_clear */
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
    const char * error_name = "strings.Error";

// splicer begin C_init_locals
// splicer end C_init_locals


    /* Create the module and add the functions */
#ifdef IS_PY3K
    m = PyModule_Create(&moduledef);
#else
    m = Py_InitModule4("strings", PY_methods,
                       PY__doc__,
                       (PyObject*)NULL,PYTHON_API_VERSION);
#endif
    if (m == NULL)
        return RETVAL;
    struct module_state *st = GETSTATE(m);


    PY_error_obj = PyErr_NewException((char *) error_name, NULL, NULL);
    if (PY_error_obj == NULL)
        return RETVAL;
    st->error = PY_error_obj;
    PyModule_AddObject(m, "Error", st->error);

// splicer begin C_init_body
// splicer end C_init_body

    /* Check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module strings");
    return RETVAL;
}
#ifdef __cplusplus
}
#endif

