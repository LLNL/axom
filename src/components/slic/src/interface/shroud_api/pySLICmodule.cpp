// pySLICmodule.cpp
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
#include "pySLICmodule.hpp"
// splicer begin include
// splicer end include

namespace asctoolkit {
namespace slic {
// splicer begin C_definition
// splicer end C_definition
PyObject *PY_error_obj;
// splicer begin additional_functions
// splicer end additional_functions

static char PY_initialize__doc__[] =
"documentation"
;

static PyObject *
PY_initialize(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.initialize
    initialize();
    Py_RETURN_NONE;
// splicer end function.initialize
}

static char PY_is_initialized__doc__[] =
"documentation"
;

static PyObject *
PY_is_initialized(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.is_initialized
    bool rv = isInitialized();
    return PyBool_FromLong(rv);
// splicer end function.is_initialized
}

static char PY_finalize__doc__[] =
"documentation"
;

static PyObject *
PY_finalize(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.finalize
    finalize();
    Py_RETURN_NONE;
// splicer end function.finalize
}
static PyMethodDef PY_methods[] = {
{"initialize", (PyCFunction)PY_initialize, METH_NOARGS, PY_initialize__doc__},
{"isInitialized", (PyCFunction)PY_is_initialized, METH_NOARGS, PY_is_initialized__doc__},
{"finalize", (PyCFunction)PY_finalize, METH_NOARGS, PY_finalize__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

/*
 * initslic - Initialization function for the module
 * *must* be called initslic
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
static int slic_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int slic_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "slic", /* m_name */
    PY__doc__, /* m_doc */
    sizeof(struct module_state), /* m_size */
    PY_methods, /* m_methods */
    NULL, /* m_reload */
    slic_traverse, /* m_traverse */
    slic_clear, /* m_clear */
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
    const char * error_name = "slic.Error";

// splicer begin C_init_locals
// splicer end C_init_locals


    /* Create the module and add the functions */
#ifdef IS_PY3K
    m = PyModule_Create(&moduledef);
#else
    m = Py_InitModule4("slic", PY_methods,
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
        Py_FatalError("can't initialize module slic");
    return RETVAL;
}
#ifdef __cplusplus
}
#endif

}  // namespace asctoolkit
}  // namespace slic
