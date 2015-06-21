// pyDataGrouptype.cpp
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
// splicer begin class.DataGroup.include
// splicer end class.DataGroup.include
namespace asctoolkit {
namespace sidre {
// splicer begin class.DataGroup.C_definition
// splicer end class.DataGroup.C_definition
// splicer begin class.DataGroup.additional_methods
// splicer end class.DataGroup.additional_methods
static int
PY_DataGroup_tp_init (PY_DataGroup *self, PyObject *args, PyObject *kwds)
{
// splicer begin class.DataGroup.type.init
    PyObject *grpobj;

    /* By requiring a PyCapsule, it is difficult to call directly from Python.
     * But the C++ constructors are private so that makes sense.
     */
    if (!PyArg_ParseTuple(args, "O!:DataGroup_init",
                          &PyCapsule_Type, &grpobj))
        return -1;

    /* capsule_dbnode */
    DataGroup *grp = static_cast<DataGroup *>(PyCapsule_GetPointer(grpobj, datagroup_capsule_name));
    self->grp = grp;
    if (grp == NULL && PyErr_Occurred())
	return -1;

    return 0;
// splicer end class.DataGroup.type.init
}

static char PY_datagroup_get_name__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_name(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.getName
const std::string & name = self->grp->getName();
PyObject * rv = PyString_FromString(name.c_str());
return rv;
// splicer end class.DataGroup.method.getName
}

static char PY_datagroup_get_parent__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_parent(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.getParent
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.getParent
}

static char PY_datagroup_get_data_store__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_data_store(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.getDataStore
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.getDataStore
}

static char PY_datagroup_get_num_views__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_num_views(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.getNumViews
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.getNumViews
}

static char PY_datagroup_get_num_groups__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_num_groups(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.getNumGroups
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.getNumGroups
}

static char PY_datagroup_has_view__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_has_view(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.hasView
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.hasView
}

static char PY_datagroup_create_view_and_buffer_simple__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_create_view_and_buffer_simple(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.createViewAndBuffer
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.createViewAndBuffer
}

static char PY_datagroup_create_view_and_buffer_from_type__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_create_view_and_buffer_from_type(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.createViewAndBuffer
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.createViewAndBuffer
}

static char PY_datagroup_create_opaque_view__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_create_opaque_view(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.createOpaqueView
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.createOpaqueView
}

static char PY_datagroup_create_view__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_create_view(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.createView
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.createView
}

static char PY_datagroup_create_external_view__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_create_external_view(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.createExternalView
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.createExternalView
}

static char PY_datagroup_move_view__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_move_view(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.moveView
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.moveView
}

static char PY_datagroup_copy_view__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_copy_view(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.copyView
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.copyView
}

static char PY_datagroup_destroy_view_and_buffer__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_destroy_view_and_buffer(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.destroyViewAndBuffer
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.destroyViewAndBuffer
}

static char PY_datagroup_get_view__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_view(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.getView
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.getView
}

static char PY_datagroup_get_view_index__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_view_index(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.getViewIndex
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.getViewIndex
}

static char PY_datagroup_get_view_name__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_view_name(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.getViewName
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.getViewName
}

static char PY_datagroup_has_group__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_has_group(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.hasGroup
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.hasGroup
}

static char PY_datagroup_create_group__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_create_group(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.createGroup
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.createGroup
}

static char PY_datagroup_move_group__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_move_group(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.moveGroup
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.moveGroup
}

static char PY_datagroup_destroy_group__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_destroy_group(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.destroyGroup
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.destroyGroup
}

static char PY_datagroup_get_group__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_group(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.getGroup
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.getGroup
}

static char PY_datagroup_get_group_index__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_group_index(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.getGroupIndex
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.getGroupIndex
}

static char PY_datagroup_get_group_name__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_group_name(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.getGroupName
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.getGroupName
}

static char PY_datagroup_print__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_print(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.print
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.print
}

static char PY_datagroup_save__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_save(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.save
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.save
}

static char PY_datagroup_load__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_load(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.load
PyErr_SetString(PyExc_NotImplementedError, "XXX");
return NULL;
// splicer end class.DataGroup.method.load
}
static PyMethodDef PY_DataGroup_methods[] = {
{"getName", (PyCFunction)PY_datagroup_get_name, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_name__doc__},
{"getParent", (PyCFunction)PY_datagroup_get_parent, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_parent__doc__},
{"getDataStore", (PyCFunction)PY_datagroup_get_data_store, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_data_store__doc__},
{"getNumViews", (PyCFunction)PY_datagroup_get_num_views, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_num_views__doc__},
{"getNumGroups", (PyCFunction)PY_datagroup_get_num_groups, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_num_groups__doc__},
{"hasView", (PyCFunction)PY_datagroup_has_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_has_view__doc__},
{"createViewAndBuffer_simple", (PyCFunction)PY_datagroup_create_view_and_buffer_simple, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_view_and_buffer_simple__doc__},
{"createViewAndBuffer_from_type", (PyCFunction)PY_datagroup_create_view_and_buffer_from_type, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_view_and_buffer_from_type__doc__},
{"createOpaqueView", (PyCFunction)PY_datagroup_create_opaque_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_opaque_view__doc__},
{"createView", (PyCFunction)PY_datagroup_create_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_view__doc__},
{"createExternalView", (PyCFunction)PY_datagroup_create_external_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_external_view__doc__},
{"moveView", (PyCFunction)PY_datagroup_move_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_move_view__doc__},
{"copyView", (PyCFunction)PY_datagroup_copy_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_copy_view__doc__},
{"destroyViewAndBuffer", (PyCFunction)PY_datagroup_destroy_view_and_buffer, METH_VARARGS|METH_KEYWORDS, PY_datagroup_destroy_view_and_buffer__doc__},
{"getView", (PyCFunction)PY_datagroup_get_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_view__doc__},
{"getViewIndex", (PyCFunction)PY_datagroup_get_view_index, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_view_index__doc__},
{"getViewName", (PyCFunction)PY_datagroup_get_view_name, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_view_name__doc__},
{"hasGroup", (PyCFunction)PY_datagroup_has_group, METH_VARARGS|METH_KEYWORDS, PY_datagroup_has_group__doc__},
{"createGroup", (PyCFunction)PY_datagroup_create_group, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_group__doc__},
{"moveGroup", (PyCFunction)PY_datagroup_move_group, METH_VARARGS|METH_KEYWORDS, PY_datagroup_move_group__doc__},
{"destroyGroup", (PyCFunction)PY_datagroup_destroy_group, METH_VARARGS|METH_KEYWORDS, PY_datagroup_destroy_group__doc__},
{"getGroup", (PyCFunction)PY_datagroup_get_group, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_group__doc__},
{"getGroupIndex", (PyCFunction)PY_datagroup_get_group_index, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_group_index__doc__},
{"getGroupName", (PyCFunction)PY_datagroup_get_group_name, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_group_name__doc__},
{"print", (PyCFunction)PY_datagroup_print, METH_VARARGS|METH_KEYWORDS, PY_datagroup_print__doc__},
{"save", (PyCFunction)PY_datagroup_save, METH_VARARGS|METH_KEYWORDS, PY_datagroup_save__doc__},
{"load", (PyCFunction)PY_datagroup_load, METH_VARARGS|METH_KEYWORDS, PY_datagroup_load__doc__},
{NULL,   (PyCFunction)NULL, 0, NULL}            /* sentinel */
};

static char DataGroup__doc__[] =
"virtual class"
;

/* static */
PyTypeObject PY_DataGroup_Type = {
        PyVarObject_HEAD_INIT(NULL, 0)
        "sidre.DataGroup",                       /* tp_name */
        sizeof(PY_DataGroup),         /* tp_basicsize */
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
        (getattrofunc)0,                /* tp_getattro */
        (setattrofunc)0,                /* tp_setattro */
        /* Functions to access object as input/output buffer */
        0,                              /* tp_as_buffer */
        /* Flags to define presence of optional/expanded features */
        Py_TPFLAGS_DEFAULT,             /* tp_flags */
        DataGroup__doc__,         /* tp_doc */
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
        PY_DataGroup_methods,                             /* tp_methods */
        0,                              /* tp_members */
        0,                             /* tp_getset */
        0,                              /* tp_base */
        0,                              /* tp_dict */
        (descrgetfunc)0,                /* tp_descr_get */
        (descrsetfunc)0,                /* tp_descr_set */
        0,                              /* tp_dictoffset */
        (initproc)PY_DataGroup_tp_init,                   /* tp_init */
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

}  // namespace asctoolkit
}  // namespace sidre
