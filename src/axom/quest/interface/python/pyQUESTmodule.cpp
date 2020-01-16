// pyQUESTmodule.cpp
// This is generated code, do not edit
//
// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "pyQUESTmodule.hpp"
#include "axom/quest/interface/inout.hpp"
#include "axom/quest/interface/signed_distance.hpp"

// splicer begin include
// splicer end include

#ifdef __cplusplus
#define SHROUD_UNUSED(param)
#else
#define SHROUD_UNUSED(param) param
#endif

#if PY_MAJOR_VERSION >= 3
#define PyInt_FromLong PyLong_FromLong
#define PyString_FromString PyUnicode_FromString
#define PyString_FromStringAndSize PyUnicode_FromStringAndSize
#endif

// splicer begin C_definition
// splicer end C_definition
PyObject* PY_error_obj;
// splicer begin additional_functions
// splicer end additional_functions

static PyObject*
PY_inout_init_mpi(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.inout_init_mpi
  const char* fileName;
  MPI_Fint comm;
  const char* SHT_kwlist[] = {
    "fileName",
    "comm",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "sO:inout_init",
                                   const_cast<char**>(SHT_kwlist), &fileName,
                                   &comm))
    return NULL;
  const std::string SH_fileName(fileName);
  MPI_Comm SH_comm = MPI_Comm_f2c(comm);
  int SHC_rv = axom::quest::inout_init(SH_fileName, SH_comm);
  PyObject* SHTPy_rv = PyInt_FromLong(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.inout_init_mpi
}

static PyObject*
PY_inout_init_serial(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.inout_init_serial
  const char* fileName;
  const char* SHT_kwlist[] = {
    "fileName",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:inout_init",
                                   const_cast<char**>(SHT_kwlist), &fileName))
    return NULL;
  const std::string SH_fileName(fileName);
  int SHC_rv = axom::quest::inout_init(SH_fileName);
  PyObject* SHTPy_rv = PyInt_FromLong(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.inout_init_serial
}

static char PY_inout_initialized__doc__[] =
  "documentation"
;

static PyObject*
PY_inout_initialized(
  PyObject* SHROUD_UNUSED(self),
  PyObject* SHROUD_UNUSED(args),
  PyObject* SHROUD_UNUSED(kwds))
{
// splicer begin function.inout_initialized
  bool SHC_rv = axom::quest::inout_initialized();
  PyObject* SHTPy_rv = PyBool_FromLong(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.inout_initialized
}

static char PY_inout_set_verbose__doc__[] =
  "documentation"
;

static PyObject*
PY_inout_set_verbose(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.inout_set_verbose
  PyObject* SHPy_verbosity;
  const char* SHT_kwlist[] = {
    "verbosity",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!:inout_set_verbose",
                                   const_cast<char**>(SHT_kwlist), &PyBool_Type,
                                   &SHPy_verbosity))
    return NULL;
  bool verbosity = PyObject_IsTrue(SHPy_verbosity);
  int SHC_rv = axom::quest::inout_set_verbose(verbosity);
  PyObject* SHTPy_rv = PyInt_FromLong(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.inout_set_verbose
}

static char PY_inout_set_vertex_weld_threshold__doc__[] =
  "documentation"
;

static PyObject*
PY_inout_set_vertex_weld_threshold(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.inout_set_vertex_weld_threshold
  double thresh;
  const char* SHT_kwlist[] = {
    "thresh",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds,
                                   "d:inout_set_vertex_weld_threshold",
                                   const_cast<char**>(SHT_kwlist), &thresh))
    return NULL;
  int SHC_rv = axom::quest::inout_set_vertex_weld_threshold(thresh);
  PyObject* SHTPy_rv = PyInt_FromLong(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.inout_set_vertex_weld_threshold
}

static char PY_inout_evaluate_1__doc__[] =
  "documentation"
;

static PyObject*
PY_inout_evaluate_1(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.inout_evaluate
  Py_ssize_t SH_nargs = 0;
  double x;
  double y;
  double z;
  const char* SHT_kwlist[] = {
    "x",
    "y",
    "z",
    NULL
  };
  bool SHC_rv;

  if (args != NULL)
    SH_nargs += PyTuple_Size(args);
  if (kwds != NULL)
    SH_nargs += PyDict_Size(args);
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "dd|d:inout_evaluate",
                                   const_cast<char**>(SHT_kwlist), &x, &y, &z))
    return NULL;
  switch (SH_nargs)
  {
  case 2:
    SHC_rv = axom::quest::inout_evaluate(x, y);
    break;
  case 3:
    SHC_rv = axom::quest::inout_evaluate(x, y, z);
    break;
  }
  PyObject* SHTPy_rv = PyBool_FromLong(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.inout_evaluate
}

static char PY_inout_get_dimension__doc__[] =
  "documentation"
;

static PyObject*
PY_inout_get_dimension(
  PyObject* SHROUD_UNUSED(self),
  PyObject* SHROUD_UNUSED(args),
  PyObject* SHROUD_UNUSED(kwds))
{
// splicer begin function.inout_get_dimension
  int SHC_rv = axom::quest::inout_get_dimension();
  PyObject* SHTPy_rv = PyInt_FromLong(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.inout_get_dimension
}

static char PY_inout_finalize__doc__[] =
  "documentation"
;

static PyObject*
PY_inout_finalize(
  PyObject* SHROUD_UNUSED(self),
  PyObject* SHROUD_UNUSED(args),
  PyObject* SHROUD_UNUSED(kwds))
{
// splicer begin function.inout_finalize
  int SHC_rv = axom::quest::inout_finalize();
  PyObject* SHTPy_rv = PyInt_FromLong(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.inout_finalize
}

static PyObject*
PY_signed_distance_init_mpi(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.signed_distance_init_mpi
  const char* file;
  MPI_Fint comm;
  const char* SHT_kwlist[] = {
    "file",
    "comm",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "sO:signed_distance_init",
                                   const_cast<char**>(SHT_kwlist), &file,
                                   &comm))
    return NULL;
  const std::string SH_file(file);
  MPI_Comm SH_comm = MPI_Comm_f2c(comm);
  int SHC_rv = axom::quest::signed_distance_init(SH_file, SH_comm);
  PyObject* SHTPy_rv = PyInt_FromLong(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.signed_distance_init_mpi
}

static PyObject*
PY_signed_distance_init_serial(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.signed_distance_init_serial
  const char* file;
  const char* SHT_kwlist[] = {
    "file",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s:signed_distance_init",
                                   const_cast<char**>(SHT_kwlist), &file))
    return NULL;
  const std::string SH_file(file);
  int SHC_rv = axom::quest::signed_distance_init(SH_file);
  PyObject* SHTPy_rv = PyInt_FromLong(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.signed_distance_init_serial
}

static char PY_signed_distance_initialized__doc__[] =
  "documentation"
;

static PyObject*
PY_signed_distance_initialized(
  PyObject* SHROUD_UNUSED(self),
  PyObject* SHROUD_UNUSED(args),
  PyObject* SHROUD_UNUSED(kwds))
{
// splicer begin function.signed_distance_initialized
  bool SHC_rv = axom::quest::signed_distance_initialized();
  PyObject* SHTPy_rv = PyBool_FromLong(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.signed_distance_initialized
}

static char PY_signed_distance_set_dimension__doc__[] =
  "documentation"
;

static PyObject*
PY_signed_distance_set_dimension(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.signed_distance_set_dimension
  int dim;
  const char* SHT_kwlist[] = {
    "dim",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds,
                                   "i:signed_distance_set_dimension",
                                   const_cast<char**>(SHT_kwlist), &dim))
    return NULL;
  axom::quest::signed_distance_set_dimension(dim);
  Py_RETURN_NONE;
// splicer end function.signed_distance_set_dimension
}

static char PY_signed_distance_set_closed_surface__doc__[] =
  "documentation"
;

static PyObject*
PY_signed_distance_set_closed_surface(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.signed_distance_set_closed_surface
  PyObject* SHPy_status;
  const char* SHT_kwlist[] = {
    "status",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds,
                                   "O!:signed_distance_set_closed_surface",
                                   const_cast<char**>(SHT_kwlist),
                                   &PyBool_Type, &SHPy_status))
    return NULL;
  bool status = PyObject_IsTrue(SHPy_status);
  axom::quest::signed_distance_set_closed_surface(status);
  Py_RETURN_NONE;
// splicer end function.signed_distance_set_closed_surface
}

static char PY_signed_distance_set_max_levels__doc__[] =
  "documentation"
;

static PyObject*
PY_signed_distance_set_max_levels(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.signed_distance_set_max_levels
  int maxLevels;
  const char* SHT_kwlist[] = {
    "maxLevels",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds,
                                   "i:signed_distance_set_max_levels",
                                   const_cast<char**>(SHT_kwlist),
                                   &maxLevels))
    return NULL;
  axom::quest::signed_distance_set_max_levels(maxLevels);
  Py_RETURN_NONE;
// splicer end function.signed_distance_set_max_levels
}

static char PY_signed_distance_set_max_occupancy__doc__[] =
  "documentation"
;

static PyObject*
PY_signed_distance_set_max_occupancy(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.signed_distance_set_max_occupancy
  int maxOccupancy;
  const char* SHT_kwlist[] = {
    "maxOccupancy",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds,
                                   "i:signed_distance_set_max_occupancy",
                                   const_cast<char**>(SHT_kwlist),
                                   &maxOccupancy))
    return NULL;
  axom::quest::signed_distance_set_max_occupancy(maxOccupancy);
  Py_RETURN_NONE;
// splicer end function.signed_distance_set_max_occupancy
}

static char PY_signed_distance_set_verbose__doc__[] =
  "documentation"
;

static PyObject*
PY_signed_distance_set_verbose(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.signed_distance_set_verbose
  PyObject* SHPy_status;
  const char* SHT_kwlist[] = {
    "status",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!:signed_distance_set_verbose",
                                   const_cast<char**>(SHT_kwlist), &PyBool_Type,
                                   &SHPy_status))
    return NULL;
  bool status = PyObject_IsTrue(SHPy_status);
  axom::quest::signed_distance_set_verbose(status);
  Py_RETURN_NONE;
// splicer end function.signed_distance_set_verbose
}

static char PY_signed_distance_use_shared_memory__doc__[] =
  "documentation"
;

static PyObject*
PY_signed_distance_use_shared_memory(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.signed_distance_use_shared_memory
  PyObject* SHPy_status;
  const char* SHT_kwlist[] = {
    "status",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds,
                                   "O!:signed_distance_use_shared_memory",
                                   const_cast<char**>(SHT_kwlist),
                                   &PyBool_Type, &SHPy_status))
    return NULL;
  bool status = PyObject_IsTrue(SHPy_status);
  axom::quest::signed_distance_use_shared_memory(status);
  Py_RETURN_NONE;
// splicer end function.signed_distance_use_shared_memory
}

static char PY_signed_distance_evaluate__doc__[] =
  "documentation"
;

static PyObject*
PY_signed_distance_evaluate(
  PyObject* SHROUD_UNUSED(self),
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.signed_distance_evaluate
  double x;
  double y;
  double z;
  const char* SHT_kwlist[] = {
    "x",
    "y",
    "z",
    NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddd:signed_distance_evaluate",
                                   const_cast<char**>(SHT_kwlist), &x, &y, &z))
    return NULL;
  double SHC_rv = axom::quest::signed_distance_evaluate(x, y, z);
  PyObject* SHTPy_rv = PyFloat_FromDouble(SHC_rv);
  return (PyObject*) SHTPy_rv;
// splicer end function.signed_distance_evaluate
}

static char PY_signed_distance_finalize__doc__[] =
  "documentation"
;

static PyObject*
PY_signed_distance_finalize(
  PyObject* SHROUD_UNUSED(self),
  PyObject* SHROUD_UNUSED(args),
  PyObject* SHROUD_UNUSED(kwds))
{
// splicer begin function.signed_distance_finalize
  axom::quest::signed_distance_finalize();
  Py_RETURN_NONE;
// splicer end function.signed_distance_finalize
}

static char PY_inout_init__doc__[] =
  "documentation"
;

static PyObject*
PY_inout_init(
  PyObject* self,
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.inout_init
  Py_ssize_t SHT_nargs = 0;
  if (args != NULL)
    SHT_nargs += PyTuple_Size(args);
  if (kwds != NULL)
    SHT_nargs += PyDict_Size(args);
  PyObject* rvobj;
  if (SHT_nargs == 2)
  {
    rvobj = PY_inout_init_mpi(self, args, kwds);
    if (!PyErr_Occurred())
    {
      return rvobj;
    }
    else if (!PyErr_ExceptionMatches(PyExc_TypeError))
    {
      return rvobj;
    }
    PyErr_Clear();
  }
  if (SHT_nargs == 1)
  {
    rvobj = PY_inout_init_serial(self, args, kwds);
    if (!PyErr_Occurred())
    {
      return rvobj;
    }
    else if (!PyErr_ExceptionMatches(PyExc_TypeError))
    {
      return rvobj;
    }
    PyErr_Clear();
  }
  PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
  return NULL;
// splicer end function.inout_init
}

static char PY_signed_distance_init__doc__[] =
  "documentation"
;

static PyObject*
PY_signed_distance_init(
  PyObject* self,
  PyObject* args,
  PyObject* kwds)
{
// splicer begin function.signed_distance_init
  Py_ssize_t SHT_nargs = 0;
  if (args != NULL)
    SHT_nargs += PyTuple_Size(args);
  if (kwds != NULL)
    SHT_nargs += PyDict_Size(args);
  PyObject* rvobj;
  if (SHT_nargs == 2)
  {
    rvobj = PY_signed_distance_init_mpi(self, args, kwds);
    if (!PyErr_Occurred())
    {
      return rvobj;
    }
    else if (!PyErr_ExceptionMatches(PyExc_TypeError))
    {
      return rvobj;
    }
    PyErr_Clear();
  }
  if (SHT_nargs == 1)
  {
    rvobj = PY_signed_distance_init_serial(self, args, kwds);
    if (!PyErr_Occurred())
    {
      return rvobj;
    }
    else if (!PyErr_ExceptionMatches(PyExc_TypeError))
    {
      return rvobj;
    }
    PyErr_Clear();
  }
  PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
  return NULL;
// splicer end function.signed_distance_init
}
static PyMethodDef PY_methods[] = {
  {"inout_initialized", (PyCFunction)PY_inout_initialized, METH_NOARGS,
   PY_inout_initialized__doc__},
  {"inout_set_verbose", (PyCFunction)PY_inout_set_verbose,
   METH_VARARGS|METH_KEYWORDS, PY_inout_set_verbose__doc__},
  {"inout_set_vertex_weld_threshold",
   (PyCFunction)PY_inout_set_vertex_weld_threshold, METH_VARARGS|METH_KEYWORDS,
   PY_inout_set_vertex_weld_threshold__doc__},
  {"inout_evaluate", (PyCFunction)PY_inout_evaluate_1,
   METH_VARARGS|METH_KEYWORDS, PY_inout_evaluate_1__doc__},
  {"inout_get_dimension", (PyCFunction)PY_inout_get_dimension, METH_NOARGS,
   PY_inout_get_dimension__doc__},
  {"inout_finalize", (PyCFunction)PY_inout_finalize, METH_NOARGS,
   PY_inout_finalize__doc__},
  {"signed_distance_initialized", (PyCFunction)PY_signed_distance_initialized,
   METH_NOARGS, PY_signed_distance_initialized__doc__},
  {"signed_distance_set_dimension",
   (PyCFunction)PY_signed_distance_set_dimension, METH_VARARGS|METH_KEYWORDS,
   PY_signed_distance_set_dimension__doc__},
  {"signed_distance_set_closed_surface",
   (PyCFunction)PY_signed_distance_set_closed_surface,
   METH_VARARGS|METH_KEYWORDS, PY_signed_distance_set_closed_surface__doc__},
  {"signed_distance_set_max_levels",
   (PyCFunction)PY_signed_distance_set_max_levels, METH_VARARGS|METH_KEYWORDS,
   PY_signed_distance_set_max_levels__doc__},
  {"signed_distance_set_max_occupancy",
   (PyCFunction)PY_signed_distance_set_max_occupancy,
   METH_VARARGS|METH_KEYWORDS, PY_signed_distance_set_max_occupancy__doc__},
  {"signed_distance_set_verbose", (PyCFunction)PY_signed_distance_set_verbose,
   METH_VARARGS|METH_KEYWORDS, PY_signed_distance_set_verbose__doc__},
  {"signed_distance_use_shared_memory",
   (PyCFunction)PY_signed_distance_use_shared_memory,
   METH_VARARGS|METH_KEYWORDS, PY_signed_distance_use_shared_memory__doc__},
  {"signed_distance_evaluate", (PyCFunction)PY_signed_distance_evaluate,
   METH_VARARGS|METH_KEYWORDS, PY_signed_distance_evaluate__doc__},
  {"signed_distance_finalize", (PyCFunction)PY_signed_distance_finalize,
   METH_NOARGS, PY_signed_distance_finalize__doc__},
  {"inout_init", (PyCFunction)PY_inout_init, METH_VARARGS|METH_KEYWORDS,
   PY_inout_init__doc__},
  {"signed_distance_init", (PyCFunction)PY_signed_distance_init,
   METH_VARARGS|METH_KEYWORDS, PY_signed_distance_init__doc__},
  {NULL,   (PyCFunction)NULL, 0, NULL}          /* sentinel */
};

/*
 * initquest - Initialization function for the module
 * *must* be called initquest
 */
static char PY__doc__[] =
  "library documentation"
;

struct module_state
{
  PyObject* error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

#if PY_MAJOR_VERSION >= 3
static int quest_traverse(PyObject* m, visitproc visit, void* arg) {
  Py_VISIT(GETSTATE(m)->error);
  return 0;
}

static int quest_clear(PyObject* m) {
  Py_CLEAR(GETSTATE(m)->error);
  return 0;
}

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "quest",   /* m_name */
  PY__doc__,   /* m_doc */
  sizeof(struct module_state),   /* m_size */
  PY_methods,   /* m_methods */
  NULL,   /* m_reload */
  quest_traverse,   /* m_traverse */
  quest_clear,   /* m_clear */
  NULL    /* m_free */
};

#define RETVAL m
#define INITERROR return NULL
#else
#define RETVAL
#define INITERROR return
#endif

extern "C" PyMODINIT_FUNC
#if PY_MAJOR_VERSION >= 3
PyInit_quest(void)
#else
initquest(void)
#endif
{
  PyObject* m = NULL;
  const char* error_name = "quest.Error";

  // splicer begin C_init_locals
  // splicer end C_init_locals


  /* Create the module and add the functions */
#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
#else
  m = Py_InitModule4("quest", PY_methods,
                     PY__doc__,
                     (PyObject*)NULL,PYTHON_API_VERSION);
#endif
  if (m == NULL)
    return RETVAL;
  struct module_state* st = GETSTATE(m);


  PY_error_obj = PyErr_NewException((char*) error_name, NULL, NULL);
  if (PY_error_obj == NULL)
    return RETVAL;
  st->error = PY_error_obj;
  PyModule_AddObject(m, "Error", st->error);

  // splicer begin C_init_body
  // splicer end C_init_body

  /* Check for errors */
  if (PyErr_Occurred())
    Py_FatalError("can't initialize module quest");
  return RETVAL;
}
