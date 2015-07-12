// pySidremodule.cpp
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
#include "pySidremodule.hpp"
// splicer begin include
// splicer end include

namespace asctoolkit {
namespace sidre {
// splicer begin C_definition
// splicer end C_definition
PyObject *PY_error_obj;
// splicer begin additional_functions
// splicer end additional_functions

static char PY_is_name_valid__doc__[] =
"documentation"
;

static PyObject *
PY_is_name_valid(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
// splicer begin function.isNameValid
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:isNameValid", kw_list,
        &name))
    {
        return NULL;
    }
    bool rv = isNameValid(name);
    return PyBool_FromLong(rv);
// splicer end function.isNameValid
}
static PyMethodDef PY_methods[] = {
{"isNameValid", (PyCFunction)PY_is_name_valid, METH_VARARGS|METH_KEYWORDS, PY_is_name_valid__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

/*
 * initsidre - Initialization function for the module
 * *must* be called initsidre
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
static int sidre_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int sidre_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "sidre", /* m_name */
    PY__doc__, /* m_doc */
    sizeof(struct module_state), /* m_size */
    PY_methods, /* m_methods */
    NULL, /* m_reload */
    sidre_traverse, /* m_traverse */
    sidre_clear, /* m_clear */
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
    const char * error_name = "sidre.Error";

// splicer begin C_init_locals
// splicer end C_init_locals


    /* Create the module and add the functions */
#ifdef IS_PY3K
    m = PyModule_Create(&moduledef);
#else
    m = Py_InitModule4("sidre", PY_methods,
                       PY__doc__,
                       (PyObject*)NULL,PYTHON_API_VERSION);
#endif
    if (m == NULL)
        return RETVAL;
    struct module_state *st = GETSTATE(m);


// DataStore
    PY_DataStore_Type.tp_new   = PyType_GenericNew;
    PY_DataStore_Type.tp_alloc = PyType_GenericAlloc;
    if (PyType_Ready(&PY_DataStore_Type) < 0)
        return RETVAL;
    Py_INCREF(&PY_DataStore_Type);
    PyModule_AddObject(m, "DataStore", (PyObject *)&PY_DataStore_Type);


// DataGroup
    PY_DataGroup_Type.tp_new   = PyType_GenericNew;
    PY_DataGroup_Type.tp_alloc = PyType_GenericAlloc;
    if (PyType_Ready(&PY_DataGroup_Type) < 0)
        return RETVAL;
    Py_INCREF(&PY_DataGroup_Type);
    PyModule_AddObject(m, "DataGroup", (PyObject *)&PY_DataGroup_Type);


// DataBuffer
    PY_DataBuffer_Type.tp_new   = PyType_GenericNew;
    PY_DataBuffer_Type.tp_alloc = PyType_GenericAlloc;
    if (PyType_Ready(&PY_DataBuffer_Type) < 0)
        return RETVAL;
    Py_INCREF(&PY_DataBuffer_Type);
    PyModule_AddObject(m, "DataBuffer", (PyObject *)&PY_DataBuffer_Type);


// DataView
    PY_DataView_Type.tp_new   = PyType_GenericNew;
    PY_DataView_Type.tp_alloc = PyType_GenericAlloc;
    if (PyType_Ready(&PY_DataView_Type) < 0)
        return RETVAL;
    Py_INCREF(&PY_DataView_Type);
    PyModule_AddObject(m, "DataView", (PyObject *)&PY_DataView_Type);


    PY_error_obj = PyErr_NewException((char *) error_name, NULL, NULL);
    if (PY_error_obj == NULL)
        return RETVAL;
    st->error = PY_error_obj;
    PyModule_AddObject(m, "Error", st->error);

// splicer begin C_init_body
// splicer end C_init_body

    /* Check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module sidre");
    return RETVAL;
}
#ifdef __cplusplus
}
#endif

}  // namespace asctoolkit
}  // namespace sidre
