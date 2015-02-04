/*
 * This is generated code.
 * Any edits must be made between the splicer.begin and splicer.end blocks.
 * All other edits will be lost.
 * Once a block is edited remove the 'UNMODIFIED' on the splicer.begin
 * comment to allow the block to be preserved when it is regenerated.
 */

#ifndef HDR_DATASTOREMODULE
#define HDR_DATASTOREMODULE
#include <Python.h>
#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif
/* DO-NOT-DELETE splicer.begin(datastore.C_declaration) MODIFIED */
#include "DatastoreInterface.hpp"
/* DO-NOT-DELETE splicer.end(datastore.C_declaration) */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/* DO-NOT-DELETE splicer.begin(datastore.class.DataObject.C_declaration) UNMODIFIED */
/* DO-NOT-DELETE splicer.end(datastore.class.DataObject.C_declaration) */

typedef struct {
    PyObject_HEAD
/* DO-NOT-DELETE splicer.begin(datastore.class.DataObject.C_object) MODIFIED */
    DataStore::DataObject *node;
/* DO-NOT-DELETE splicer.end(datastore.class.DataObject.C_object) */
} DS_DataObjectObject;

extern PyTypeObject DS_DataObject_Type;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/* DO-NOT-DELETE splicer.begin(datastore.class.Function.C_declaration) UNMODIFIED */
/* DO-NOT-DELETE splicer.end(datastore.class.Function.C_declaration) */

typedef struct {
    PyObject_HEAD
/* DO-NOT-DELETE splicer.begin(datastore.class.Function.C_object) MODIFIED */
    DataStore::DataObject *node;
/* DO-NOT-DELETE splicer.end(datastore.class.Function.C_object) */
} DS_FunctionObject;

extern PyTypeObject DS_Function_Type;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/* DO-NOT-DELETE splicer.begin(datastore.class.DataGroup.C_declaration) UNMODIFIED */
/* DO-NOT-DELETE splicer.end(datastore.class.DataGroup.C_declaration) */

typedef struct {
    PyObject_HEAD
/* DO-NOT-DELETE splicer.begin(datastore.class.DataGroup.C_object) MODIFIED */
    DataStore::DataGroup *node;
/* DO-NOT-DELETE splicer.end(datastore.class.DataGroup.C_object) */
} DS_DataGroupObject;

extern PyTypeObject DS_DataGroup_Type;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

extern PyObject *DS_error_obj;

#ifdef IS_PY3K
#define MOD_INITDATASTORE PyInit_datastore
#else
#define MOD_INITDATASTORE initdatastore
#endif
PyMODINIT_FUNC MOD_INITDATASTORE(void);
#endif
