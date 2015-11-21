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
// splicer begin class.DataGroup.impl.include
// splicer end class.DataGroup.impl.include
namespace asctoolkit {
namespace sidre {
// splicer begin class.DataGroup.impl.C_definition
// splicer end class.DataGroup.impl.C_definition
// splicer begin class.DataGroup.impl.additional_methods
// splicer end class.DataGroup.impl.additional_methods
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
    DataGroup *grp = static_cast<DataGroup *>(PyCapsule_GetPointer(grpobj, PY_DataGroup_capsule_name));
    self->BBB = grp;
    if (grp == NULL && PyErr_Occurred())
	return -1;

    return 0;
// splicer end class.DataGroup.type.init
}
static PyObject *
PY_DataGroup_tp_richcompare (PY_DataGroup *self, PyObject *other, int opid)
{
// splicer begin class.DataGroup.type.richcompare
PyObject *rv = Py_NotImplemented;
if (PyObject_IsInstance(other, (PyObject*) &PY_DataGroup_Type)) {
    PY_DataGroup *pyother = (PY_DataGroup *) other;
    switch (opid) {
    case Py_EQ:
	if (self->BBB == pyother->BBB) {
	    rv = Py_True;
	} else {
	    rv = Py_False;
	}
	break;
    case Py_NE:
	if (self->BBB != pyother->BBB) {
	    rv = Py_True;
	} else {
	    rv = Py_False;
	}
	break;
    case Py_LT:
    case Py_LE:
    case Py_GE:
    case Py_GT:
	break;
    }
 }
Py_INCREF(rv);
return rv;
// splicer end class.DataGroup.type.richcompare
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
// splicer begin class.DataGroup.method.get_name
    const std::string & rv = self->BBB->getName();
    return PyString_FromString(rv.c_str());
// splicer end class.DataGroup.method.get_name
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
// splicer begin class.DataGroup.method.get_parent
    DataGroup * rv = self->BBB->getParent();
    PY_DataGroup * rv_obj = PyObject_New(PY_DataGroup, &PY_DataGroup_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.get_parent
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
// splicer begin class.DataGroup.method.get_data_store
    DataStore * rv = self->BBB->getDataStore();
    PY_DataStore * rv_obj = PyObject_New(PY_DataStore, &PY_DataStore_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.get_data_store
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
// splicer begin class.DataGroup.method.get_num_views
    size_t rv = self->BBB->getNumViews();
    return PyInt_FromLong(rv);
// splicer end class.DataGroup.method.get_num_views
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
// splicer begin class.DataGroup.method.get_num_groups
    size_t rv = self->BBB->getNumGroups();
    return PyInt_FromLong(rv);
// splicer end class.DataGroup.method.get_num_groups
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
// splicer begin class.DataGroup.method.has_view
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:hasView", kw_list,
        &name))
    {
        return NULL;
    }
    bool rv = self->BBB->hasView(name);
    return PyBool_FromLong(rv);
// splicer end class.DataGroup.method.has_view
}

static char PY_datagroup_get_view_from_name__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_view_from_name(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.get_view_from_name
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:getView", kw_list,
        &name))
    {
        return NULL;
    }
    DataView * rv = self->BBB->getView(name);
    PY_DataView * rv_obj = PyObject_New(PY_DataView, &PY_DataView_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.get_view_from_name
}

static char PY_datagroup_get_view_from_index__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_get_view_from_index(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.get_view_from_index
    const ATK_IndexType idx;
    const char *kwcpp = "idx";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:getView", kw_list,
        &idx))
    {
        return NULL;
    }
    DataView * rv = self->BBB->getView(idx);
    PY_DataView * rv_obj = PyObject_New(PY_DataView, &PY_DataView_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.get_view_from_index
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
// splicer begin class.DataGroup.method.get_view_index
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:getViewIndex", kw_list,
        &name))
    {
        return NULL;
    }
    IndexType rv = self->BBB->getViewIndex(name);
    return Py_BuildValue("i", rv);
// splicer end class.DataGroup.method.get_view_index
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
// splicer begin class.DataGroup.method.get_view_name
    ATK_IndexType idx;
    const char *kwcpp = "idx";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:getViewName", kw_list,
        &idx))
    {
        return NULL;
    }
    const std::string & rv = self->BBB->getViewName(idx);
    if (! nameIsValid(rv)) {
        Py_RETURN_NONE;
    }
    
    return PyString_FromString(rv.c_str());
// splicer end class.DataGroup.method.get_view_name
}

static char PY_datagroup_create_view_and_allocate_from_type__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_create_view_and_allocate_from_type(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.create_view_and_allocate_from_type
    const char * name;
    int type;
    ATK_SidreLength numelems;
    const char *kwcpp = "name\0type\0numelems";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5,(char *) kwcpp+10, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "sil:createViewAndAllocate", kw_list,
        &name, &type, &numelems))
    {
        return NULL;
    }
    DataView * rv = self->BBB->createViewAndAllocate(name, getTypeID(type), numelems);
    PY_DataView * rv_obj = PyObject_New(PY_DataView, &PY_DataView_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.create_view_and_allocate_from_type
}

static char PY_datagroup_create_view_empty__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_create_view_empty(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.create_view_empty
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:createView", kw_list,
        &name))
    {
        return NULL;
    }
    DataView * rv = self->BBB->createView(name);
    PY_DataView * rv_obj = PyObject_New(PY_DataView, &PY_DataView_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.create_view_empty
}

static char PY_datagroup_create_view_into_buffer__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_create_view_into_buffer(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.create_view_into_buffer
    const char * name;
    PY_DataBuffer * buff;
    DataBuffer * buff_ptr;
    const char *kwcpp = "name\0buff";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "sO!:createView", kw_list,
        &name, &PY_DataBuffer_Type, &buff))
    {
        return NULL;
    }
    buff_ptr = (buff ? buff->BBB : NULL);
    DataView * rv = self->BBB->createView(name, buff_ptr);
    PY_DataView * rv_obj = PyObject_New(PY_DataView, &PY_DataView_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.create_view_into_buffer
}

static char PY_datagroup_create_view_into_buffer_nelems__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_create_view_into_buffer_nelems(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.create_view_into_buffer_nelems
    const char * name;
    PY_DataBuffer * buff;
    DataBuffer * buff_ptr;
    ATK_SidreLength numelems;
    const char *kwcpp = "name\0buff\0numelems";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5,(char *) kwcpp+10, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "sO!l:createView", kw_list,
        &name, &PY_DataBuffer_Type, &buff, &numelems))
    {
        return NULL;
    }
    buff_ptr = (buff ? buff->BBB : NULL);
    DataView * rv = self->BBB->createView(name, buff_ptr, numelems);
    PY_DataView * rv_obj = PyObject_New(PY_DataView, &PY_DataView_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.create_view_into_buffer_nelems
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
// splicer begin class.DataGroup.method.create_opaque_view
    const char * name;
    void * opaque_ptr;
    const char *kwcpp = "name\0opaque_ptr";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "sO:createOpaqueView", kw_list,
        &name, &opaque_ptr))
    {
        return NULL;
    }
    DataView * rv = self->BBB->createOpaqueView(name, opaque_ptr);
    PY_DataView * rv_obj = PyObject_New(PY_DataView, &PY_DataView_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.create_opaque_view
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
// splicer begin class.DataGroup.method.create_external_view
    const char * name;
    void * external_data;
    int type;
    ATK_SidreLength numelems;
    const char *kwcpp = "name\0external_data\0type\0numelems";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+5,(char *) kwcpp+19,(char *) kwcpp+24, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "sOil:createExternalView", kw_list,
        &name, &external_data, &type, &numelems))
    {
        return NULL;
    }
    DataView * rv = self->BBB->createExternalView(name, external_data, getTypeID(type), numelems);
    PY_DataView * rv_obj = PyObject_New(PY_DataView, &PY_DataView_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.create_external_view
}

static char PY_datagroup_destroy_view__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_destroy_view(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.destroy_view
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:destroyView", kw_list,
        &name))
    {
        return NULL;
    }
    self->BBB->destroyView(name);
    Py_RETURN_NONE;
// splicer end class.DataGroup.method.destroy_view
}

static char PY_datagroup_destroy_view_and_data__doc__[] =
"documentation"
;

static PyObject *
PY_datagroup_destroy_view_and_data(
  PY_DataGroup *self,
  PyObject *args,
  PyObject *kwds)
{
// splicer begin class.DataGroup.method.destroy_view_and_data
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:destroyViewAndData", kw_list,
        &name))
    {
        return NULL;
    }
    self->BBB->destroyViewAndData(name);
    Py_RETURN_NONE;
// splicer end class.DataGroup.method.destroy_view_and_data
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
// splicer begin class.DataGroup.method.move_view
    PY_DataView * view;
    DataView * view_ptr;
    const char *kwcpp = "view";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!:moveView", kw_list,
        &PY_DataView_Type, &view))
    {
        return NULL;
    }
    view_ptr = (view ? view->BBB : NULL);
    DataView * rv = self->BBB->moveView(view_ptr);
    PY_DataView * rv_obj = PyObject_New(PY_DataView, &PY_DataView_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.move_view
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
// splicer begin class.DataGroup.method.copy_view
    PY_DataView * view;
    DataView * view_ptr;
    const char *kwcpp = "view";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!:copyView", kw_list,
        &PY_DataView_Type, &view))
    {
        return NULL;
    }
    view_ptr = (view ? view->BBB : NULL);
    DataView * rv = self->BBB->copyView(view_ptr);
    PY_DataView * rv_obj = PyObject_New(PY_DataView, &PY_DataView_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.copy_view
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
// splicer begin class.DataGroup.method.has_group
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:hasGroup", kw_list,
        &name))
    {
        return NULL;
    }
    bool rv = self->BBB->hasGroup(name);
    return PyBool_FromLong(rv);
// splicer end class.DataGroup.method.has_group
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
// splicer begin class.DataGroup.method.get_group
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:getGroup", kw_list,
        &name))
    {
        return NULL;
    }
    DataGroup * rv = self->BBB->getGroup(name);
    if (rv == ATK_NULLPTR) {
        Py_RETURN_NONE;
    }
    
    PY_DataGroup * rv_obj = PyObject_New(PY_DataGroup, &PY_DataGroup_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.get_group
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
// splicer begin class.DataGroup.method.get_group_index
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:getGroupIndex", kw_list,
        &name))
    {
        return NULL;
    }
    IndexType rv = self->BBB->getGroupIndex(name);
    return Py_BuildValue("i", rv);
// splicer end class.DataGroup.method.get_group_index
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
// splicer begin class.DataGroup.method.get_group_name
    ATK_IndexType idx;
    const char *kwcpp = "idx";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "i:getGroupName", kw_list,
        &idx))
    {
        return NULL;
    }
    const std::string & rv = self->BBB->getGroupName(idx);
    if (! nameIsValid(rv)) {
        Py_RETURN_NONE;
    }
    
    return PyString_FromString(rv.c_str());
// splicer end class.DataGroup.method.get_group_name
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
// splicer begin class.DataGroup.method.create_group
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:createGroup", kw_list,
        &name))
    {
        return NULL;
    }
    DataGroup * rv = self->BBB->createGroup(name);
    PY_DataGroup * rv_obj = PyObject_New(PY_DataGroup, &PY_DataGroup_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.create_group
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
// splicer begin class.DataGroup.method.destroy_group
    const char * name;
    const char *kwcpp = "name";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:destroyGroup", kw_list,
        &name))
    {
        return NULL;
    }
    self->BBB->destroyGroup(name);
    Py_RETURN_NONE;
// splicer end class.DataGroup.method.destroy_group
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
// splicer begin class.DataGroup.method.move_group
    PY_DataGroup * grp;
    DataGroup * grp_ptr;
    const char *kwcpp = "grp";
    char *kw_list[] = { (char *) kwcpp+0, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!:moveGroup", kw_list,
        &PY_DataGroup_Type, &grp))
    {
        return NULL;
    }
    grp_ptr = (grp ? grp->BBB : NULL);
    DataGroup * rv = self->BBB->moveGroup(grp_ptr);
    PY_DataGroup * rv_obj = PyObject_New(PY_DataGroup, &PY_DataGroup_Type);
    rv_obj->BBB = rv;
    return (PyObject *) rv_obj;
// splicer end class.DataGroup.method.move_group
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
    self->BBB->print();
    Py_RETURN_NONE;
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
    const char * obase;
    const char * protocol;
    const char *kwcpp = "obase\0protocol";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+6, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ss:save", kw_list,
        &obase, &protocol))
    {
        return NULL;
    }
    self->BBB->save(obase, protocol);
    Py_RETURN_NONE;
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
    const char * obase;
    const char * protocol;
    const char *kwcpp = "obase\0protocol";
    char *kw_list[] = { (char *) kwcpp+0,(char *) kwcpp+6, NULL };
    
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ss:load", kw_list,
        &obase, &protocol))
    {
        return NULL;
    }
    self->BBB->load(obase, protocol);
    Py_RETURN_NONE;
// splicer end class.DataGroup.method.load
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
// splicer begin class.DataGroup.get_view
    int numNamedArgs = (kwds ? PyDict_Size(kwds) : 0);
    int numArgs = PyTuple_GET_SIZE(args);
    int totArgs = numArgs + numNamedArgs;
    PyObject *rvobj;
    {
        rvobj = PY_datagroup_get_view_from_name(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    {
        rvobj = PY_datagroup_get_view_from_index(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
    return NULL;
// splicer end class.DataGroup.get_view
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
// splicer begin class.DataGroup.create_view
    int numNamedArgs = (kwds ? PyDict_Size(kwds) : 0);
    int numArgs = PyTuple_GET_SIZE(args);
    int totArgs = numArgs + numNamedArgs;
    PyObject *rvobj;
    {
        rvobj = PY_datagroup_create_view_empty(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    {
        rvobj = PY_datagroup_create_view_into_buffer(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    {
        rvobj = PY_datagroup_create_view_into_buffer_nelems(self, args, kwds);
        if (!PyErr_Occurred()) {
            return rvobj;
        } else if (! PyErr_ExceptionMatches(PyExc_TypeError)) {
            return rvobj;
        }
        PyErr_Clear();
    }
    PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
    return NULL;
// splicer end class.DataGroup.create_view
}
// splicer begin class.DataGroup.impl.after_methods
// splicer end class.DataGroup.impl.after_methods
static PyMethodDef PY_DataGroup_methods[] = {
{"getName", (PyCFunction)PY_datagroup_get_name, METH_NOARGS, PY_datagroup_get_name__doc__},
{"getParent", (PyCFunction)PY_datagroup_get_parent, METH_NOARGS, PY_datagroup_get_parent__doc__},
{"getDataStore", (PyCFunction)PY_datagroup_get_data_store, METH_NOARGS, PY_datagroup_get_data_store__doc__},
{"getNumViews", (PyCFunction)PY_datagroup_get_num_views, METH_NOARGS, PY_datagroup_get_num_views__doc__},
{"getNumGroups", (PyCFunction)PY_datagroup_get_num_groups, METH_NOARGS, PY_datagroup_get_num_groups__doc__},
{"hasView", (PyCFunction)PY_datagroup_has_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_has_view__doc__},
{"getView_from_name", (PyCFunction)PY_datagroup_get_view_from_name, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_view_from_name__doc__},
{"getView_from_index", (PyCFunction)PY_datagroup_get_view_from_index, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_view_from_index__doc__},
{"getViewIndex", (PyCFunction)PY_datagroup_get_view_index, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_view_index__doc__},
{"getViewName", (PyCFunction)PY_datagroup_get_view_name, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_view_name__doc__},
{"createViewAndAllocate_from_type", (PyCFunction)PY_datagroup_create_view_and_allocate_from_type, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_view_and_allocate_from_type__doc__},
{"createView_empty", (PyCFunction)PY_datagroup_create_view_empty, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_view_empty__doc__},
{"createView_into_buffer", (PyCFunction)PY_datagroup_create_view_into_buffer, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_view_into_buffer__doc__},
{"createView_into_buffer_nelems", (PyCFunction)PY_datagroup_create_view_into_buffer_nelems, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_view_into_buffer_nelems__doc__},
{"createOpaqueView", (PyCFunction)PY_datagroup_create_opaque_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_opaque_view__doc__},
{"createExternalView", (PyCFunction)PY_datagroup_create_external_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_external_view__doc__},
{"destroyView", (PyCFunction)PY_datagroup_destroy_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_destroy_view__doc__},
{"destroyViewAndData", (PyCFunction)PY_datagroup_destroy_view_and_data, METH_VARARGS|METH_KEYWORDS, PY_datagroup_destroy_view_and_data__doc__},
{"moveView", (PyCFunction)PY_datagroup_move_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_move_view__doc__},
{"copyView", (PyCFunction)PY_datagroup_copy_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_copy_view__doc__},
{"hasGroup", (PyCFunction)PY_datagroup_has_group, METH_VARARGS|METH_KEYWORDS, PY_datagroup_has_group__doc__},
{"getGroup", (PyCFunction)PY_datagroup_get_group, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_group__doc__},
{"getGroupIndex", (PyCFunction)PY_datagroup_get_group_index, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_group_index__doc__},
{"getGroupName", (PyCFunction)PY_datagroup_get_group_name, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_group_name__doc__},
{"createGroup", (PyCFunction)PY_datagroup_create_group, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_group__doc__},
{"destroyGroup", (PyCFunction)PY_datagroup_destroy_group, METH_VARARGS|METH_KEYWORDS, PY_datagroup_destroy_group__doc__},
{"moveGroup", (PyCFunction)PY_datagroup_move_group, METH_VARARGS|METH_KEYWORDS, PY_datagroup_move_group__doc__},
{"print", (PyCFunction)PY_datagroup_print, METH_NOARGS, PY_datagroup_print__doc__},
{"save", (PyCFunction)PY_datagroup_save, METH_VARARGS|METH_KEYWORDS, PY_datagroup_save__doc__},
{"load", (PyCFunction)PY_datagroup_load, METH_VARARGS|METH_KEYWORDS, PY_datagroup_load__doc__},
{"getView", (PyCFunction)PY_datagroup_get_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_get_view__doc__},
{"createView", (PyCFunction)PY_datagroup_create_view, METH_VARARGS|METH_KEYWORDS, PY_datagroup_create_view__doc__},
// splicer begin class.DataGroup.PyMethodDef
// splicer end class.DataGroup.PyMethodDef
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
        (richcmpfunc)PY_DataGroup_tp_richcompare,                 /* tp_richcompare */
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
