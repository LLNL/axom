// pyQUESTmodule.cpp
// This file is generated by Shroud 0.13.0. Do not edit.
//
// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "pyQUESTmodule.hpp"

// splicer begin include
// splicer end include

#ifdef __cplusplus
  #define SHROUD_UNUSED(param)
#else
  #define SHROUD_UNUSED(param) param
#endif

#if PY_MAJOR_VERSION >= 3
  #define PyInt_AsLong PyLong_AsLong
  #define PyInt_FromLong PyLong_FromLong
  #define PyInt_FromSize_t PyLong_FromSize_t
  #define PyString_FromString PyUnicode_FromString
  #define PyString_FromStringAndSize PyUnicode_FromStringAndSize
#endif

// splicer begin C_definition
// splicer end C_definition
PyObject *PY_error_obj;
// splicer begin additional_functions
// splicer end additional_functions

#ifdef AXOM_USE_MPI
static PyObject *PY_inout_init_mpi(PyObject *SHROUD_UNUSED(self),
                                   PyObject *args,
                                   PyObject *kwds)
{
  // splicer begin function.inout_init_mpi
  char *fileName;
  MPI_Fint comm;
  const char *SHT_kwlist[] = {"fileName", "comm", nullptr};
  PyObject *SHTPy_rv = nullptr;

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "sO:inout_init",
                                  const_cast<char **>(SHT_kwlist),
                                  &fileName,
                                  &comm))
    return nullptr;
  const std::string SH_fileName(fileName);
  MPI_Comm SH_comm = MPI_Comm_f2c(comm);
  int SHCXX_rv = axom::quest::inout_init(SH_fileName, SH_comm);
  SHTPy_rv = PyInt_FromLong(SHCXX_rv);
  return (PyObject *)SHTPy_rv;
  // splicer end function.inout_init_mpi
}
#endif  // ifdef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
static PyObject *PY_inout_init_serial(PyObject *SHROUD_UNUSED(self),
                                      PyObject *args,
                                      PyObject *kwds)
{
  // splicer begin function.inout_init_serial
  char *fileName;
  const char *SHT_kwlist[] = {"fileName", nullptr};
  PyObject *SHTPy_rv = nullptr;

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "s:inout_init",
                                  const_cast<char **>(SHT_kwlist),
                                  &fileName))
    return nullptr;
  const std::string SH_fileName(fileName);
  int SHCXX_rv = axom::quest::inout_init(SH_fileName);
  SHTPy_rv = PyInt_FromLong(SHCXX_rv);
  return (PyObject *)SHTPy_rv;
  // splicer end function.inout_init_serial
}
#endif  // ifndef AXOM_USE_MPI

static char PY_inout_initialized__doc__[] = "documentation";

static PyObject *PY_inout_initialized(PyObject *SHROUD_UNUSED(self),
                                      PyObject *SHROUD_UNUSED(args),
                                      PyObject *SHROUD_UNUSED(kwds))
{
  // splicer begin function.inout_initialized
  PyObject *SHTPy_rv = nullptr;

  bool SHCXX_rv = axom::quest::inout_initialized();
  SHTPy_rv = PyBool_FromLong(SHCXX_rv);
  if(SHTPy_rv == nullptr) goto fail;
  return (PyObject *)SHTPy_rv;

fail:
  Py_XDECREF(SHTPy_rv);
  return nullptr;
  // splicer end function.inout_initialized
}

static char PY_inout_set_dimension__doc__[] = "documentation";

static PyObject *PY_inout_set_dimension(PyObject *SHROUD_UNUSED(self),
                                        PyObject *args,
                                        PyObject *kwds)
{
  // splicer begin function.inout_set_dimension
  int dim;
  const char *SHT_kwlist[] = {"dim", nullptr};
  PyObject *SHTPy_rv = nullptr;

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "i:inout_set_dimension",
                                  const_cast<char **>(SHT_kwlist),
                                  &dim))
    return nullptr;
  int SHCXX_rv = axom::quest::inout_set_dimension(dim);
  SHTPy_rv = PyInt_FromLong(SHCXX_rv);
  return (PyObject *)SHTPy_rv;
  // splicer end function.inout_set_dimension
}

static char PY_inout_set_verbose__doc__[] = "documentation";

static PyObject *PY_inout_set_verbose(PyObject *SHROUD_UNUSED(self),
                                      PyObject *args,
                                      PyObject *kwds)
{
  // splicer begin function.inout_set_verbose
  bool verbosity;
  PyObject *SHPy_verbosity;
  const char *SHT_kwlist[] = {"verbosity", nullptr};
  PyObject *SHTPy_rv = nullptr;

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "O!:inout_set_verbose",
                                  const_cast<char **>(SHT_kwlist),
                                  &PyBool_Type,
                                  &SHPy_verbosity))
    return nullptr;
  verbosity = PyObject_IsTrue(SHPy_verbosity);
  int SHCXX_rv = axom::quest::inout_set_verbose(verbosity);
  SHTPy_rv = PyInt_FromLong(SHCXX_rv);
  return (PyObject *)SHTPy_rv;
  // splicer end function.inout_set_verbose
}

static char PY_inout_set_vertex_weld_threshold__doc__[] = "documentation";

static PyObject *PY_inout_set_vertex_weld_threshold(PyObject *SHROUD_UNUSED(self),
                                                    PyObject *args,
                                                    PyObject *kwds)
{
  // splicer begin function.inout_set_vertex_weld_threshold
  double thresh;
  const char *SHT_kwlist[] = {"thresh", nullptr};
  PyObject *SHTPy_rv = nullptr;

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "d:inout_set_vertex_weld_threshold",
                                  const_cast<char **>(SHT_kwlist),
                                  &thresh))
    return nullptr;
  int SHCXX_rv = axom::quest::inout_set_vertex_weld_threshold(thresh);
  SHTPy_rv = PyInt_FromLong(SHCXX_rv);
  return (PyObject *)SHTPy_rv;
  // splicer end function.inout_set_vertex_weld_threshold
}

static char PY_inout_set_segments_per_knot_span__doc__[] = "documentation";

static PyObject *PY_inout_set_segments_per_knot_span(PyObject *SHROUD_UNUSED(self),
                                                     PyObject *args,
                                                     PyObject *kwds)
{
  // splicer begin function.inout_set_segments_per_knot_span
  int segmentsPerKnotSpan;
  const char *SHT_kwlist[] = {"segmentsPerKnotSpan", nullptr};
  PyObject *SHTPy_rv = nullptr;

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "i:inout_set_segments_per_knot_span",
                                  const_cast<char **>(SHT_kwlist),
                                  &segmentsPerKnotSpan))
    return nullptr;
  int SHCXX_rv =
    axom::quest::inout_set_segments_per_knot_span(segmentsPerKnotSpan);
  SHTPy_rv = PyInt_FromLong(SHCXX_rv);
  return (PyObject *)SHTPy_rv;
  // splicer end function.inout_set_segments_per_knot_span
}

static char PY_inout_evaluate_1__doc__[] = "documentation";

static PyObject *PY_inout_evaluate_1(PyObject *SHROUD_UNUSED(self),
                                     PyObject *args,
                                     PyObject *kwds)
{
  // splicer begin function.inout_evaluate
  Py_ssize_t SH_nargs = 0;
  double x;
  double y;
  double z;
  const char *SHT_kwlist[] = {"x", "y", "z", nullptr};
  bool SHCXX_rv;
  PyObject *SHTPy_rv = nullptr;

  if(args != nullptr) SH_nargs += PyTuple_Size(args);
  if(kwds != nullptr) SH_nargs += PyDict_Size(args);
  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "dd|d:inout_evaluate",
                                  const_cast<char **>(SHT_kwlist),
                                  &x,
                                  &y,
                                  &z))
    return nullptr;
  switch(SH_nargs)
  {
  case 2:
    SHCXX_rv = axom::quest::inout_evaluate(x, y);
    break;
  case 3:
    SHCXX_rv = axom::quest::inout_evaluate(x, y, z);
    break;
  default:
    PyErr_SetString(PyExc_ValueError, "Wrong number of arguments");
    return nullptr;
  }
  SHTPy_rv = PyBool_FromLong(SHCXX_rv);
  if(SHTPy_rv == nullptr) goto fail;
  return (PyObject *)SHTPy_rv;

fail:
  Py_XDECREF(SHTPy_rv);
  return nullptr;
  // splicer end function.inout_evaluate
}

static char PY_inout_get_dimension__doc__[] = "documentation";

static PyObject *PY_inout_get_dimension(PyObject *SHROUD_UNUSED(self),
                                        PyObject *SHROUD_UNUSED(args),
                                        PyObject *SHROUD_UNUSED(kwds))
{
  // splicer begin function.inout_get_dimension
  PyObject *SHTPy_rv = nullptr;

  int SHCXX_rv = axom::quest::inout_get_dimension();
  SHTPy_rv = PyInt_FromLong(SHCXX_rv);
  return (PyObject *)SHTPy_rv;
  // splicer end function.inout_get_dimension
}

static char PY_inout_finalize__doc__[] = "documentation";

static PyObject *PY_inout_finalize(PyObject *SHROUD_UNUSED(self),
                                   PyObject *SHROUD_UNUSED(args),
                                   PyObject *SHROUD_UNUSED(kwds))
{
  // splicer begin function.inout_finalize
  PyObject *SHTPy_rv = nullptr;

  int SHCXX_rv = axom::quest::inout_finalize();
  SHTPy_rv = PyInt_FromLong(SHCXX_rv);
  return (PyObject *)SHTPy_rv;
  // splicer end function.inout_finalize
}

#ifdef AXOM_USE_MPI
static PyObject *PY_signed_distance_init_mpi(PyObject *SHROUD_UNUSED(self),
                                             PyObject *args,
                                             PyObject *kwds)
{
  // splicer begin function.signed_distance_init_mpi
  char *file;
  MPI_Fint comm;
  const char *SHT_kwlist[] = {"file", "comm", nullptr};
  PyObject *SHTPy_rv = nullptr;

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "sO:signed_distance_init",
                                  const_cast<char **>(SHT_kwlist),
                                  &file,
                                  &comm))
    return nullptr;
  const std::string SH_file(file);
  MPI_Comm SH_comm = MPI_Comm_f2c(comm);
  int SHCXX_rv = axom::quest::signed_distance_init(SH_file, SH_comm);
  SHTPy_rv = PyInt_FromLong(SHCXX_rv);
  return (PyObject *)SHTPy_rv;
  // splicer end function.signed_distance_init_mpi
}
#endif  // ifdef AXOM_USE_MPI

#ifndef AXOM_USE_MPI
static PyObject *PY_signed_distance_init_serial(PyObject *SHROUD_UNUSED(self),
                                                PyObject *args,
                                                PyObject *kwds)
{
  // splicer begin function.signed_distance_init_serial
  char *file;
  const char *SHT_kwlist[] = {"file", nullptr};
  PyObject *SHTPy_rv = nullptr;

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "s:signed_distance_init",
                                  const_cast<char **>(SHT_kwlist),
                                  &file))
    return nullptr;
  const std::string SH_file(file);
  int SHCXX_rv = axom::quest::signed_distance_init(SH_file);
  SHTPy_rv = PyInt_FromLong(SHCXX_rv);
  return (PyObject *)SHTPy_rv;
  // splicer end function.signed_distance_init_serial
}
#endif  // ifndef AXOM_USE_MPI

static char PY_signed_distance_initialized__doc__[] = "documentation";

static PyObject *PY_signed_distance_initialized(PyObject *SHROUD_UNUSED(self),
                                                PyObject *SHROUD_UNUSED(args),
                                                PyObject *SHROUD_UNUSED(kwds))
{
  // splicer begin function.signed_distance_initialized
  PyObject *SHTPy_rv = nullptr;

  bool SHCXX_rv = axom::quest::signed_distance_initialized();
  SHTPy_rv = PyBool_FromLong(SHCXX_rv);
  if(SHTPy_rv == nullptr) goto fail;
  return (PyObject *)SHTPy_rv;

fail:
  Py_XDECREF(SHTPy_rv);
  return nullptr;
  // splicer end function.signed_distance_initialized
}

static char PY_signed_distance_set_dimension__doc__[] = "documentation";

static PyObject *PY_signed_distance_set_dimension(PyObject *SHROUD_UNUSED(self),
                                                  PyObject *args,
                                                  PyObject *kwds)
{
  // splicer begin function.signed_distance_set_dimension
  int dim;
  const char *SHT_kwlist[] = {"dim", nullptr};

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "i:signed_distance_set_dimension",
                                  const_cast<char **>(SHT_kwlist),
                                  &dim))
    return nullptr;
  axom::quest::signed_distance_set_dimension(dim);
  Py_RETURN_NONE;
  // splicer end function.signed_distance_set_dimension
}

static char PY_signed_distance_set_closed_surface__doc__[] = "documentation";

static PyObject *PY_signed_distance_set_closed_surface(PyObject *SHROUD_UNUSED(self),
                                                       PyObject *args,
                                                       PyObject *kwds)
{
  // splicer begin function.signed_distance_set_closed_surface
  bool status;
  PyObject *SHPy_status;
  const char *SHT_kwlist[] = {"status", nullptr};

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "O!:signed_distance_set_closed_surface",
                                  const_cast<char **>(SHT_kwlist),
                                  &PyBool_Type,
                                  &SHPy_status))
    return nullptr;
  status = PyObject_IsTrue(SHPy_status);
  axom::quest::signed_distance_set_closed_surface(status);
  Py_RETURN_NONE;
  // splicer end function.signed_distance_set_closed_surface
}

static char PY_signed_distance_set_compute_signs__doc__[] = "documentation";

static PyObject *PY_signed_distance_set_compute_signs(PyObject *SHROUD_UNUSED(self),
                                                      PyObject *args,
                                                      PyObject *kwds)
{
  // splicer begin function.signed_distance_set_compute_signs
  bool computeSign;
  PyObject *SHPy_computeSign;
  const char *SHT_kwlist[] = {"computeSign", nullptr};

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "O!:signed_distance_set_compute_signs",
                                  const_cast<char **>(SHT_kwlist),
                                  &PyBool_Type,
                                  &SHPy_computeSign))
    return nullptr;
  computeSign = PyObject_IsTrue(SHPy_computeSign);
  axom::quest::signed_distance_set_compute_signs(computeSign);
  Py_RETURN_NONE;
  // splicer end function.signed_distance_set_compute_signs
}

static char PY_signed_distance_set_allocator__doc__[] = "documentation";

static PyObject *PY_signed_distance_set_allocator(PyObject *SHROUD_UNUSED(self),
                                                  PyObject *args,
                                                  PyObject *kwds)
{
  // splicer begin function.signed_distance_set_allocator
  int allocatorID;
  const char *SHT_kwlist[] = {"allocatorID", nullptr};

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "i:signed_distance_set_allocator",
                                  const_cast<char **>(SHT_kwlist),
                                  &allocatorID))
    return nullptr;
  axom::quest::signed_distance_set_allocator(allocatorID);
  Py_RETURN_NONE;
  // splicer end function.signed_distance_set_allocator
}

static char PY_signed_distance_set_verbose__doc__[] = "documentation";

static PyObject *PY_signed_distance_set_verbose(PyObject *SHROUD_UNUSED(self),
                                                PyObject *args,
                                                PyObject *kwds)
{
  // splicer begin function.signed_distance_set_verbose
  bool status;
  PyObject *SHPy_status;
  const char *SHT_kwlist[] = {"status", nullptr};

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "O!:signed_distance_set_verbose",
                                  const_cast<char **>(SHT_kwlist),
                                  &PyBool_Type,
                                  &SHPy_status))
    return nullptr;
  status = PyObject_IsTrue(SHPy_status);
  axom::quest::signed_distance_set_verbose(status);
  Py_RETURN_NONE;
  // splicer end function.signed_distance_set_verbose
}

static char PY_signed_distance_use_shared_memory__doc__[] = "documentation";

static PyObject *PY_signed_distance_use_shared_memory(PyObject *SHROUD_UNUSED(self),
                                                      PyObject *args,
                                                      PyObject *kwds)
{
  // splicer begin function.signed_distance_use_shared_memory
  bool status;
  PyObject *SHPy_status;
  const char *SHT_kwlist[] = {"status", nullptr};

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "O!:signed_distance_use_shared_memory",
                                  const_cast<char **>(SHT_kwlist),
                                  &PyBool_Type,
                                  &SHPy_status))
    return nullptr;
  status = PyObject_IsTrue(SHPy_status);
  axom::quest::signed_distance_use_shared_memory(status);
  Py_RETURN_NONE;
  // splicer end function.signed_distance_use_shared_memory
}

static char PY_signed_distance_set_execution_space__doc__[] = "documentation";

static PyObject *PY_signed_distance_set_execution_space(PyObject *SHROUD_UNUSED(self),
                                                        PyObject *args,
                                                        PyObject *kwds)
{
  // splicer begin function.signed_distance_set_execution_space
  int execSpace;
  const char *SHT_kwlist[] = {"execSpace", nullptr};

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "i:signed_distance_set_execution_space",
                                  const_cast<char **>(SHT_kwlist),
                                  &execSpace))
    return nullptr;
  axom::quest::SignedDistExec SH_execSpace =
    static_cast<axom::quest::SignedDistExec>(execSpace);
  axom::quest::signed_distance_set_execution_space(SH_execSpace);
  Py_RETURN_NONE;
  // splicer end function.signed_distance_set_execution_space
}

static PyObject *PY_signed_distance_evaluate_0(PyObject *SHROUD_UNUSED(self),
                                               PyObject *args,
                                               PyObject *kwds)
{
  // splicer begin function.signed_distance_evaluate_0
  double x;
  double y;
  double z;
  const char *SHT_kwlist[] = {"x", "y", "z", nullptr};
  PyObject *SHTPy_rv = nullptr;

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "ddd:signed_distance_evaluate",
                                  const_cast<char **>(SHT_kwlist),
                                  &x,
                                  &y,
                                  &z))
    return nullptr;
  double SHCXX_rv = axom::quest::signed_distance_evaluate(x, y, z);
  SHTPy_rv = PyFloat_FromDouble(SHCXX_rv);
  return (PyObject *)SHTPy_rv;
  // splicer end function.signed_distance_evaluate_0
}

static PyObject *PY_signed_distance_evaluate_1(PyObject *SHROUD_UNUSED(self),
                                               PyObject *args,
                                               PyObject *kwds)
{
  // splicer begin function.signed_distance_evaluate_1
  double x;
  double y;
  double z;
  double cp_x;
  double cp_y;
  double cp_z;
  double n_x;
  double n_y;
  double n_z;
  const char *SHT_kwlist[] =
    {"x", "y", "z", "cp_x", "cp_y", "cp_z", "n_x", "n_y", "n_z", nullptr};
  PyObject *SHTPy_rv = nullptr;  // return value object

  if(!PyArg_ParseTupleAndKeywords(args,
                                  kwds,
                                  "ddddddddd:signed_distance_evaluate",
                                  const_cast<char **>(SHT_kwlist),
                                  &x,
                                  &y,
                                  &z,
                                  &cp_x,
                                  &cp_y,
                                  &cp_z,
                                  &n_x,
                                  &n_y,
                                  &n_z))
    return nullptr;
  double SHCXX_rv =
    axom::quest::signed_distance_evaluate(x, y, z, cp_x, cp_y, cp_z, n_x, n_y, n_z);
  SHTPy_rv = Py_BuildValue("ddddddd", SHCXX_rv, cp_x, cp_y, cp_z, n_x, n_y, n_z);
  return SHTPy_rv;
  // splicer end function.signed_distance_evaluate_1
}

static char PY_signed_distance_finalize__doc__[] = "documentation";

static PyObject *PY_signed_distance_finalize(PyObject *SHROUD_UNUSED(self),
                                             PyObject *SHROUD_UNUSED(args),
                                             PyObject *SHROUD_UNUSED(kwds))
{
  // splicer begin function.signed_distance_finalize
  axom::quest::signed_distance_finalize();
  Py_RETURN_NONE;
  // splicer end function.signed_distance_finalize
}

static char PY_inout_init__doc__[] = "documentation";

static PyObject *PY_inout_init(PyObject *self, PyObject *args, PyObject *kwds)
{
  // splicer begin function.inout_init
  Py_ssize_t SHT_nargs = 0;
  if(args != nullptr) SHT_nargs += PyTuple_Size(args);
  if(kwds != nullptr) SHT_nargs += PyDict_Size(args);
  PyObject *rvobj;
#ifdef AXOM_USE_MPI
  if(SHT_nargs == 2)
  {
    rvobj = PY_inout_init_mpi(self, args, kwds);
    if(!PyErr_Occurred())
    {
      return rvobj;
    }
    else if(!PyErr_ExceptionMatches(PyExc_TypeError))
    {
      return rvobj;
    }
    PyErr_Clear();
  }
#endif  // ifdef AXOM_USE_MPI
#ifndef AXOM_USE_MPI
  if(SHT_nargs == 1)
  {
    rvobj = PY_inout_init_serial(self, args, kwds);
    if(!PyErr_Occurred())
    {
      return rvobj;
    }
    else if(!PyErr_ExceptionMatches(PyExc_TypeError))
    {
      return rvobj;
    }
    PyErr_Clear();
  }
#endif  // ifndef AXOM_USE_MPI
  PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
  return nullptr;
  // splicer end function.inout_init
}

static char PY_signed_distance_init__doc__[] = "documentation";

static PyObject *PY_signed_distance_init(PyObject *self,
                                         PyObject *args,
                                         PyObject *kwds)
{
  // splicer begin function.signed_distance_init
  Py_ssize_t SHT_nargs = 0;
  if(args != nullptr) SHT_nargs += PyTuple_Size(args);
  if(kwds != nullptr) SHT_nargs += PyDict_Size(args);
  PyObject *rvobj;
#ifdef AXOM_USE_MPI
  if(SHT_nargs == 2)
  {
    rvobj = PY_signed_distance_init_mpi(self, args, kwds);
    if(!PyErr_Occurred())
    {
      return rvobj;
    }
    else if(!PyErr_ExceptionMatches(PyExc_TypeError))
    {
      return rvobj;
    }
    PyErr_Clear();
  }
#endif  // ifdef AXOM_USE_MPI
#ifndef AXOM_USE_MPI
  if(SHT_nargs == 1)
  {
    rvobj = PY_signed_distance_init_serial(self, args, kwds);
    if(!PyErr_Occurred())
    {
      return rvobj;
    }
    else if(!PyErr_ExceptionMatches(PyExc_TypeError))
    {
      return rvobj;
    }
    PyErr_Clear();
  }
#endif  // ifndef AXOM_USE_MPI
  PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
  return nullptr;
  // splicer end function.signed_distance_init
}

static char PY_signed_distance_evaluate__doc__[] = "documentation";

static PyObject *PY_signed_distance_evaluate(PyObject *self,
                                             PyObject *args,
                                             PyObject *kwds)
{
  // splicer begin function.signed_distance_evaluate
  Py_ssize_t SHT_nargs = 0;
  if(args != nullptr) SHT_nargs += PyTuple_Size(args);
  if(kwds != nullptr) SHT_nargs += PyDict_Size(args);
  PyObject *rvobj;
  if(SHT_nargs == 3)
  {
    rvobj = PY_signed_distance_evaluate_0(self, args, kwds);
    if(!PyErr_Occurred())
    {
      return rvobj;
    }
    else if(!PyErr_ExceptionMatches(PyExc_TypeError))
    {
      return rvobj;
    }
    PyErr_Clear();
  }
  if(SHT_nargs == 9)
  {
    rvobj = PY_signed_distance_evaluate_1(self, args, kwds);
    if(!PyErr_Occurred())
    {
      return rvobj;
    }
    else if(!PyErr_ExceptionMatches(PyExc_TypeError))
    {
      return rvobj;
    }
    PyErr_Clear();
  }
  PyErr_SetString(PyExc_TypeError, "wrong arguments multi-dispatch");
  return nullptr;
  // splicer end function.signed_distance_evaluate
}
static PyMethodDef PY_methods[] = {
  {"inout_initialized",
   (PyCFunction)PY_inout_initialized,
   METH_NOARGS,
   PY_inout_initialized__doc__},
  {"inout_set_dimension",
   (PyCFunction)PY_inout_set_dimension,
   METH_VARARGS | METH_KEYWORDS,
   PY_inout_set_dimension__doc__},
  {"inout_set_verbose",
   (PyCFunction)PY_inout_set_verbose,
   METH_VARARGS | METH_KEYWORDS,
   PY_inout_set_verbose__doc__},
  {"inout_set_vertex_weld_threshold",
   (PyCFunction)PY_inout_set_vertex_weld_threshold,
   METH_VARARGS | METH_KEYWORDS,
   PY_inout_set_vertex_weld_threshold__doc__},
  {"inout_set_segments_per_knot_span",
   (PyCFunction)PY_inout_set_segments_per_knot_span,
   METH_VARARGS | METH_KEYWORDS,
   PY_inout_set_segments_per_knot_span__doc__},
  {"inout_evaluate",
   (PyCFunction)PY_inout_evaluate_1,
   METH_VARARGS | METH_KEYWORDS,
   PY_inout_evaluate_1__doc__},
  {"inout_get_dimension",
   (PyCFunction)PY_inout_get_dimension,
   METH_NOARGS,
   PY_inout_get_dimension__doc__},
  {"inout_finalize",
   (PyCFunction)PY_inout_finalize,
   METH_NOARGS,
   PY_inout_finalize__doc__},
  {"signed_distance_initialized",
   (PyCFunction)PY_signed_distance_initialized,
   METH_NOARGS,
   PY_signed_distance_initialized__doc__},
  {"signed_distance_set_dimension",
   (PyCFunction)PY_signed_distance_set_dimension,
   METH_VARARGS | METH_KEYWORDS,
   PY_signed_distance_set_dimension__doc__},
  {"signed_distance_set_closed_surface",
   (PyCFunction)PY_signed_distance_set_closed_surface,
   METH_VARARGS | METH_KEYWORDS,
   PY_signed_distance_set_closed_surface__doc__},
  {"signed_distance_set_compute_signs",
   (PyCFunction)PY_signed_distance_set_compute_signs,
   METH_VARARGS | METH_KEYWORDS,
   PY_signed_distance_set_compute_signs__doc__},
  {"signed_distance_set_allocator",
   (PyCFunction)PY_signed_distance_set_allocator,
   METH_VARARGS | METH_KEYWORDS,
   PY_signed_distance_set_allocator__doc__},
  {"signed_distance_set_verbose",
   (PyCFunction)PY_signed_distance_set_verbose,
   METH_VARARGS | METH_KEYWORDS,
   PY_signed_distance_set_verbose__doc__},
  {"signed_distance_use_shared_memory",
   (PyCFunction)PY_signed_distance_use_shared_memory,
   METH_VARARGS | METH_KEYWORDS,
   PY_signed_distance_use_shared_memory__doc__},
  {"signed_distance_set_execution_space",
   (PyCFunction)PY_signed_distance_set_execution_space,
   METH_VARARGS | METH_KEYWORDS,
   PY_signed_distance_set_execution_space__doc__},
  {"signed_distance_finalize",
   (PyCFunction)PY_signed_distance_finalize,
   METH_NOARGS,
   PY_signed_distance_finalize__doc__},
  {"inout_init",
   (PyCFunction)PY_inout_init,
   METH_VARARGS | METH_KEYWORDS,
   PY_inout_init__doc__},
  {"signed_distance_init",
   (PyCFunction)PY_signed_distance_init,
   METH_VARARGS | METH_KEYWORDS,
   PY_signed_distance_init__doc__},
  {"signed_distance_evaluate",
   (PyCFunction)PY_signed_distance_evaluate,
   METH_VARARGS | METH_KEYWORDS,
   PY_signed_distance_evaluate__doc__},
  {nullptr, (PyCFunction) nullptr, 0, nullptr} /* sentinel */
};

/*
 * initquest - Initialization function for the module
 * *must* be called initquest
 */
static char PY__doc__[] = "library documentation";

struct module_state
{
  PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
  #define GETSTATE(m) ((struct module_state *)PyModule_GetState(m))
#else
  #define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

#if PY_MAJOR_VERSION >= 3
static int quest_traverse(PyObject *m, visitproc visit, void *arg)
{
  Py_VISIT(GETSTATE(m)->error);
  return 0;
}

static int quest_clear(PyObject *m)
{
  Py_CLEAR(GETSTATE(m)->error);
  return 0;
}

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "quest",                     /* m_name */
  PY__doc__,                   /* m_doc */
  sizeof(struct module_state), /* m_size */
  PY_methods,                  /* m_methods */
  nullptr,                     /* m_reload */
  quest_traverse,              /* m_traverse */
  quest_clear,                 /* m_clear */
  NULL                         /* m_free */
};

  #define RETVAL m
  #define INITERROR return nullptr
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
  PyObject *m = nullptr;
  const char *error_name = "quest.Error";

  // splicer begin C_init_locals
  // splicer end C_init_locals

  /* Create the module and add the functions */
#if PY_MAJOR_VERSION >= 3
  m = PyModule_Create(&moduledef);
#else
  m = Py_InitModule4("quest",
                     PY_methods,
                     PY__doc__,
                     (PyObject *)nullptr,
                     PYTHON_API_VERSION);
#endif
  if(m == nullptr) return RETVAL;
  struct module_state *st = GETSTATE(m);

  // enum axom::quest::SignedDistExec
  PyModule_AddIntConstant(m, "CPU", axom::quest::CPU);
  PyModule_AddIntConstant(m, "OpenMP", axom::quest::OpenMP);
  PyModule_AddIntConstant(m, "GPU", axom::quest::GPU);

  PY_error_obj = PyErr_NewException((char *)error_name, nullptr, nullptr);
  if(PY_error_obj == nullptr) return RETVAL;
  st->error = PY_error_obj;
  PyModule_AddObject(m, "Error", st->error);

  // splicer begin C_init_body
  // splicer end C_init_body

  /* Check for errors */
  if(PyErr_Occurred()) Py_FatalError("can't initialize module quest");
  return RETVAL;
}
