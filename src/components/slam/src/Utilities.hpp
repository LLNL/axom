/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/**
 * \file
 * \brief A few utility functions used by the MeshAPI component.
 */
#ifndef MESHAPI_UTILITIES_H_
#define MESHAPI_UTILITIES_H_

#include <string>
#include <cmath>

namespace asctoolkit {
namespace meshapi {

  typedef int MeshIndexType;
  typedef int MeshSizeType;

  class NotImplementedException {};

namespace util {


/** \brief A helper class to print the name of a few types */
  template<typename T> struct TypeToString { static std::string to_string(){return "<unspecialized>"; }
  };

/** \brief A helper class to print the name of integers as 'int' */
  template<> struct TypeToString<int>{ static std::string       to_string(){return "int"; }
  };

/** \brief A helper class to print the name of doubles as 'double' */
  template<> struct TypeToString<double>{ static std::string    to_string(){return "double"; }
  };


} // end namespace util
} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_UTILITIES_H_
