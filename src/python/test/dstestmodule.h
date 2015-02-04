/*
 * This is generated code.
 * Any edits must be made between the splicer.begin and splicer.end blocks.
 * All other edits will be lost.
 * Once a block is edited remove the 'UNMODIFIED' on the splicer.begin
 * comment to allow the block to be preserved when it is regenerated.
 */

#ifndef HDR_DSTESTMODULE
#define HDR_DSTESTMODULE
#include <Python.h>
#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif
/* DO-NOT-DELETE splicer.begin(dstest.C_declaration) UNMODIFIED */
/* DO-NOT-DELETE splicer.end(dstest.C_declaration) */

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

extern PyObject *DSTST_error_obj;

#ifdef IS_PY3K
#define MOD_INITDSTEST PyInit_dstest
#else
#define MOD_INITDSTEST initdstest
#endif
PyMODINIT_FUNC MOD_INITDSTEST(void);
#endif
