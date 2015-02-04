/*
 * This is generated code.
 * Any edits must be made between the splicer.begin and splicer.end blocks.
 * All other edits will be lost.
 * Once a block is edited remove the 'UNMODIFIED' on the splicer.begin
 * comment to allow the block to be preserved when it is regenerated.
 */

#include <dsmodule.h>
/* DO-NOT-DELETE splicer.begin(datastore.C_definition) MODIFIED */
#include "DataFunctionPython.hpp"

static PyObject *capsule_arg(void *voidptr)
{
    PyObject *voidobj = PyCapsule_New(voidptr, "DS_capsule", NULL);
    PyObject *cargs = PyTuple_New(1);
    PyTuple_SET_ITEM(cargs, 0, voidobj);
    return cargs;
}

/* DO-NOT-DELETE splicer.end(datastore.C_definition) */
PyObject *DS_error_obj;

/*--------------------------------------------------------------------------*/

static char DS_DataObject_name__doc__[] = "";

static PyObject *
DS_DataObject_name_get(DS_DataObjectObject *self, void *context)
{
/* DO-NOT-DELETE splicer.begin(datastore.class.DataObject.descriptor.name.get) UNMODIFIED */
    PyErr_SetString(PyExc_NotImplementedError, "DS_DataObject_name_get");
    return NULL;
/* DO-NOT-DELETE splicer.end(datastore.class.DataObject.descriptor.name.get) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static PyGetSetDef DS_DataObject_getset[] = {
{"name", (getter) DS_DataObject_name_get, NULL, DS_DataObject_name__doc__, NULL},
{NULL}            /* sentinel */
};

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static PyObject *
DS_DataObject_repr(DS_DataObjectObject *self)
{
/* DO-NOT-DELETE splicer.begin(datastore.class.DataObject.type.repr) MODIFIED */
    char buffer[100];  // XXX - buffer size
    snprintf(buffer, 100, "<%s %s>",
             Py_TYPE(self)->tp_name, self->node->Name().c_str());
    buffer[99] = '\0';
    return PyUnicode_FromString(buffer);
/* DO-NOT-DELETE splicer.end(datastore.class.DataObject.type.repr) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static int
DS_DataObject_init(DS_DataObjectObject *self, PyObject *args, PyObject *kwds)
{
/* DO-NOT-DELETE splicer.begin(datastore.class.DataObject.type.init) MODIFIED */
    PyObject *dsobj;

    /* By requiring a PyCapsule, it is difficult to call directly
     * from Python.
     * Intended to be called from compiled code via PyObject_Call
     */
    if (!PyArg_ParseTuple(args, "O!:DS_DataObject_init",
                          &PyCapsule_Type, &dsobj))
        return -1;

    auto node = static_cast<DataStore::DataObject *>(PyCapsule_GetPointer(dsobj, "DS_capsule"));
    self->node = node;
    if (node == nullptr && PyErr_Occurred())
	return -1;

    return 0;
/* DO-NOT-DELETE splicer.end(datastore.class.DataObject.type.init) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static char DS_DataObject__doc__[] =
"base object"
;

/* static */
PyTypeObject DS_DataObject_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "datastore.DataObject",                       /* tp_name */
        sizeof(DS_DataObjectObject),         /* tp_basicsize */
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
        (reprfunc)DS_DataObject_repr,                    /* tp_repr */
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
        DS_DataObject__doc__,         /* tp_doc */
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
        0,                             /* tp_methods */
        0,                              /* tp_members */
        DS_DataObject_getset,                             /* tp_getset */
        0,                              /* tp_base */
        0,                              /* tp_dict */
        (descrgetfunc)0,                /* tp_descr_get */
        (descrsetfunc)0,                /* tp_descr_set */
        0,                              /* tp_dictoffset */
        (initproc)DS_DataObject_init,                   /* tp_init */
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

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static PyObject *
DS_Function_call(DS_FunctionObject *self, PyObject *args, PyObject *kwds)
{
/* DO-NOT-DELETE splicer.begin(datastore.class.Function.type.call) MODIFIED */
    auto *func = self->node->GetData<DataStore::DataFunction*>();

    // XXX deal with arguments  PyObject to DataArray?
    func->Call();

    Py_RETURN_NONE;
/* DO-NOT-DELETE splicer.end(datastore.class.Function.type.call) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static char DS_Function__doc__[] =
"datastore function"
;

/* static */
PyTypeObject DS_Function_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "datastore.Function",                       /* tp_name */
        sizeof(DS_FunctionObject),         /* tp_basicsize */
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
        (ternaryfunc)DS_Function_call,                 /* tp_call */
        (reprfunc)0,                    /* tp_str */
        (getattrofunc)0,                /* tp_getattro */
        (setattrofunc)0,                /* tp_setattro */
        /* Functions to access object as input/output buffer */
        0,                              /* tp_as_buffer */
        /* Flags to define presence of optional/expanded features */
        Py_TPFLAGS_DEFAULT,             /* tp_flags */
        DS_Function__doc__,         /* tp_doc */
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
        0,                             /* tp_methods */
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

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static char DS_DataGroup___dir____doc__[] = "";

static PyObject *
DS_DataGroup___dir__(
  DS_DataGroupObject *self,
  PyObject *args,
  PyObject *kwds)
{
/* DO-NOT-DELETE splicer.begin(datastore.class.DataGroup.method.__dir__) UNMODIFIED */
    PyErr_SetString(PyExc_NotImplementedError, "DS_DataGroup___dir__");
    return NULL;
/* DO-NOT-DELETE splicer.end(datastore.class.DataGroup.method.__dir__) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static PyMethodDef DS_DataGroup_methods[] = {
{"__dir__", (PyCFunction)DS_DataGroup___dir__, METH_VARARGS|METH_KEYWORDS, DS_DataGroup___dir____doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static char DS_DataGroup___doc____doc__[] = "";


/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static PyGetSetDef DS_DataGroup_getset[] = {
{"__doc__", NULL, NULL, DS_DataGroup___doc____doc__, NULL},
{NULL}            /* sentinel */
};

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static PyObject *
DS_DataGroup_getattro(DS_DataGroupObject *self, PyObject *name)
{
/* DO-NOT-DELETE splicer.begin(datastore.class.DataGroup.type.getattro) MODIFIED */
    PyObject *value = nullptr;

    if (PyUnicode_Check(name)) {
        const char *sname = PyUnicode_AsUTF8(name);
	DataStore::DataObject* const node = self->node->GetDataObject(sname);
        if (node != NULL) {
	    PyObject *cargs = capsule_arg(node);
	    value = PyObject_Call((PyObject *) &DS_DataObject_Type, cargs, NULL);
	    Py_DECREF(cargs);
        }
    }
    
    if (value == nullptr) {
        /* Fall back on generic to get __class__ and __dict__ */
        return PyObject_GenericGetAttr((PyObject *)self, name);
    }

    return value;
/* DO-NOT-DELETE splicer.end(datastore.class.DataGroup.type.getattro) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static int
DS_DataGroup_setattro(DS_DataGroupObject *self, PyObject *name, PyObject *value)
{
/* DO-NOT-DELETE splicer.begin(datastore.class.DataGroup.type.setattro) UNMODIFIED */
    PyErr_SetString(PyExc_NotImplementedError, "setattro");
    return -1;
/* DO-NOT-DELETE splicer.end(datastore.class.DataGroup.type.setattro) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static char DS_DataGroup__doc__[] =
"DataStore Group"
;

/* static */
PyTypeObject DS_DataGroup_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "datastore.DataGroup",                       /* tp_name */
        sizeof(DS_DataGroupObject),         /* tp_basicsize */
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
        (getattrofunc)DS_DataGroup_getattro,                /* tp_getattro */
        (setattrofunc)DS_DataGroup_setattro,                /* tp_setattro */
        /* Functions to access object as input/output buffer */
        0,                              /* tp_as_buffer */
        /* Flags to define presence of optional/expanded features */
        Py_TPFLAGS_DEFAULT,             /* tp_flags */
        DS_DataGroup__doc__,         /* tp_doc */
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
        DS_DataGroup_methods,                             /* tp_methods */
        0,                              /* tp_members */
        DS_DataGroup_getset,                             /* tp_getset */
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

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static char DS_create_data_store__doc__[] =
"Create a datastore"
;

static PyObject *
DS_create_data_store(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
/* DO-NOT-DELETE splicer.begin(datastore.function.create_data_store) MODIFIED */
    char *kw_list[] = { "name", NULL };
    const char *name;

    if (!PyArg_ParseTupleAndKeywords(
            args, kwds, "s:create_data_store", kw_list,
            &name))
        return NULL;

    DataStore::DataGroup* const grp = DataStore::CreateDataStore(name);
    assert(grp != NULL);

    PyObject *cargs = capsule_arg(grp);
    PyObject *rv = PyObject_Call((PyObject *) &DS_DataGroup_Type, cargs, NULL);
    Py_DECREF(cargs);

    return rv;
    //    Py_RETURN_NONE;
/* DO-NOT-DELETE splicer.end(datastore.function.create_data_store) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static char DS_add_object__doc__[] =
"Add object to a group."
;

static PyObject *
DS_add_object(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
/* DO-NOT-DELETE splicer.begin(datastore.function.add_object) MODIFIED */
    PyObject *dataobj;
    DS_DataGroupObject *dsgrp;
    DataStore::DataGroup *grp;
    const char *name;

    if (!PyArg_ParseTuple(args, "O!sO:add_object",
                          &DS_DataGroup_Type, &dsgrp,
			  &name, &dataobj))
	return NULL;

    auto *node = dsgrp->node->CreateDataObject(name);

    PyObject *cargs = capsule_arg(node);
    PyObject *rv = PyObject_Call((PyObject *) &DS_DataObject_Type, cargs, NULL);
    Py_DECREF(cargs);

    return rv;
/* DO-NOT-DELETE splicer.end(datastore.function.add_object) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static char DS_add_function__doc__[] =
"Add function to a group."
;

static PyObject *
DS_add_function(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
/* DO-NOT-DELETE splicer.begin(datastore.function.add_function) MODIFIED */
    PyObject *funcobj;
    DS_DataGroupObject *dsgrp;
    DataStore::DataGroup *grp;
    const char *name;

    if (!PyArg_ParseTuple(args, "O!sO:add_object",
                          &DS_DataGroup_Type, &dsgrp,
			  &name, &funcobj))
	return NULL;

    auto *func = new DataStore::DataFunctionPython(funcobj);
    auto *node = dsgrp->node->CreateDataObject(name)
	->SetDataPointer(func);

    PyObject *cargs = capsule_arg(node);
    PyObject *rv = PyObject_Call((PyObject *) &DS_Function_Type, cargs, NULL);
    Py_DECREF(cargs);

    return rv;
/* DO-NOT-DELETE splicer.end(datastore.function.add_function) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static char DS_stop_here__doc__[] =
"A function to help debug"
;

static PyObject *
DS_stop_here(
  PyObject *self,    /* not used */
  PyObject *args,
  PyObject *kwds)
{
/* DO-NOT-DELETE splicer.begin(datastore.function.stop_here) UNMODIFIED */
    PyErr_SetString(PyExc_NotImplementedError, "DS_stop_here");
    return NULL;
/* DO-NOT-DELETE splicer.end(datastore.function.stop_here) */
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

static PyMethodDef DS_methods[] = {
{"create_data_store", (PyCFunction)DS_create_data_store, METH_VARARGS|METH_KEYWORDS, DS_create_data_store__doc__},
{"add_object", (PyCFunction)DS_add_object, METH_VARARGS|METH_KEYWORDS, DS_add_object__doc__},
{"add_function", (PyCFunction)DS_add_function, METH_VARARGS|METH_KEYWORDS, DS_add_function__doc__},
{"stop_here", (PyCFunction)DS_stop_here, METH_VARARGS|METH_KEYWORDS, DS_stop_here__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*
 * initdatastore - Initialization function for the module
 * *must* be called initdatastore
 */
static char DS__doc__[] =
"Datastore"
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
static int datastore_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int datastore_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "datastore", /* m_name */
    DS__doc__, /* m_doc */
    sizeof(struct module_state), /* m_size */
    DS_methods, /* m_methods */
    NULL, /* m_reload */
    datastore_traverse, /* m_traverse */
    datastore_clear, /* m_clear */
    NULL  /* m_free */
};

#define RETVAL m
#define INITERROR return NULL

PyMODINIT_FUNC
PyInit_datastore(void)

#else
#define RETVAL
#define INITERROR return

PyMODINIT_FUNC
initdatastore(void)
#endif
{
    PyObject *m = NULL;
/* DO-NOT-DELETE splicer.begin(datastore.C_init_locals) UNMODIFIED */
/* DO-NOT-DELETE splicer.end(datastore.C_init_locals) */

    /* Create the module and add the functions */
#ifdef IS_PY3K
    m = PyModule_Create(&moduledef);
#else
    m = Py_InitModule4("datastore", DS_methods,
                       DS__doc__,
                       (PyObject*)NULL,PYTHON_API_VERSION);
#endif
    if (m == NULL)
        return RETVAL;
    struct module_state *st = GETSTATE(m);

    DS_DataObject_Type.tp_new   = PyType_GenericNew;
    DS_DataObject_Type.tp_alloc = PyType_GenericAlloc;
    if (PyType_Ready(&DS_DataObject_Type) < 0)
        return RETVAL;
    Py_INCREF(&DS_DataObject_Type);
    PyModule_AddObject(m, "DataObject", (PyObject *)&DS_DataObject_Type);

    DS_Function_Type.tp_new   = PyType_GenericNew;
    DS_Function_Type.tp_alloc = PyType_GenericAlloc;
    DS_Function_Type.tp_base  = &DS_DataObject_Type;
    if (PyType_Ready(&DS_Function_Type) < 0)
        return RETVAL;
    Py_INCREF(&DS_Function_Type);
    PyModule_AddObject(m, "Function", (PyObject *)&DS_Function_Type);

    DS_DataGroup_Type.tp_new   = PyType_GenericNew;
    DS_DataGroup_Type.tp_alloc = PyType_GenericAlloc;
    DS_DataGroup_Type.tp_base  = &DS_DataObject_Type;
    if (PyType_Ready(&DS_DataGroup_Type) < 0)
        return RETVAL;
    Py_INCREF(&DS_DataGroup_Type);
    PyModule_AddObject(m, "DataGroup", (PyObject *)&DS_DataGroup_Type);

    DS_error_obj = PyErr_NewException("datastore.Error", NULL, NULL);
    if (DS_error_obj == NULL)
        return RETVAL;
    st->error = DS_error_obj;
    PyModule_AddObject(m, "Error", st->error);

/* DO-NOT-DELETE splicer.begin(datastore.C_init_body) UNMODIFIED */
/* DO-NOT-DELETE splicer.end(datastore.C_init_body) */

    /* Check for errors */
    if (PyErr_Occurred())
        Py_FatalError("can't initialize module datastore");
    return RETVAL;
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

