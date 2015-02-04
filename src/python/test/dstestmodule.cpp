/*
 * This is generated code.
 * Any edits must be made between the splicer.begin and splicer.end blocks.
 * All other edits will be lost.
 * Once a block is edited remove the 'UNMODIFIED' on the splicer.begin
 * comment to allow the block to be preserved when it is regenerated.
 */

#include <dstestmodule.h>
/* DO-NOT-DELETE splicer.begin(dstest.C_definition) MODIFIED */
#include "dsmodule.h"
#include "DatastoreInterface.hpp"
#include "DataFunction.hpp"
/* DO-NOT-DELETE splicer.end(dstest.C_definition) */
PyObject *DSTST_error_obj;

/*--------------------------------------------------------------------------*/

static char DSTST_call_python_function__doc__[] =
"call a Python funtion via C"
;

static PyObject *
DSTST_call_python_function(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
/* DO-NOT-DELETE splicer.begin(dstest.function.call_python_function) MODIFIED */
    /*
     * Pass in DataObject, extract DataFunction, call
     */
    DS_FunctionObject *dataobj;
    char *kw_list[] = { "function", NULL };
    
    const char *name;

    if (!PyArg_ParseTupleAndKeywords(
            args, kwds, "O!:call_python_function", kw_list,
	    &DS_Function_Type, &dataobj))
        return NULL;

    auto *func = dataobj->node->GetData<DataStore::DataFunction*>();
    func->Call();

    Py_RETURN_NONE;
/* DO-NOT-DELETE splicer.end(dstest.function.call_python_function) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static PyMethodDef DSTST_methods[] = {
{"call_python_function", (PyCFunction)DSTST_call_python_function, METH_VARARGS|METH_KEYWORDS, DSTST_call_python_function__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*
 * initdstest - Initialization function for the module
 * *must* be called initdstest
 */
static char DSTST__doc__[] =
"Datastore test functions"
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
static int dstest_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int dstest_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "dstest", /* m_name */
    DSTST__doc__, /* m_doc */
    sizeof(struct module_state), /* m_size */
    DSTST_methods, /* m_methods */
    NULL, /* m_reload */
    dstest_traverse, /* m_traverse */
    dstest_clear, /* m_clear */
    NULL  /* m_free */
};

#define RETVAL m
#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_dstest(void)

#else
#define RETVAL
#define INITERROR return

PyMODINIT_FUNC
initdstest(void)
#endif
{
    PyObject *m = NULL;
/* DO-NOT-DELETE splicer.begin(dstest.C_init_locals) UNMODIFIED */
/* DO-NOT-DELETE splicer.end(dstest.C_init_locals) */

    /* Create the module and add the functions */
#ifdef IS_PY3K
    m = PyModule_Create(&moduledef);
#else
    m = Py_InitModule4("dstest", DSTST_methods,
                       DSTST__doc__,
                       (PyObject*)NULL,PYTHON_API_VERSION);
#endif
    if (m == NULL)
        return RETVAL;
    struct module_state *st = GETSTATE(m);

    DSTST_error_obj = PyErr_NewException("dstest.Error", NULL, NULL);
    if (DSTST_error_obj == NULL)
        return RETVAL;
    st->error = DSTST_error_obj;
    PyModule_AddObject(m, "Error", st->error);

/* DO-NOT-DELETE splicer.begin(dstest.C_init_body) UNMODIFIED */
/* DO-NOT-DELETE splicer.end(dstest.C_init_body) */

    /* Check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module dstest");
    return RETVAL;
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

