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
 * \brief A few utility functions used by the SLAM component.
 */
#ifndef SLAM_UTILITIES_H_
#define SLAM_UTILITIES_H_

#include <string>

namespace asctoolkit {
namespace slam {

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


  /**
   * \brief Simple utility function to search through ancestor directories to find a file
   * \param [in] fileName The path to the original file.
   * \param [in] numAncestors The maximum number of ancestor directories to try
   * \return A string representing a valid path to the file when this can be found,
   *         else, will return the original fileName.
   */
  std::string findFileInAncestorDirs(const std::string& fileName, int numAncestors = 4);

} // end namespace util
} // end namespace slam
} // end namespace asctoolkit

#endif //  SLAM_UTILITIES_H_
