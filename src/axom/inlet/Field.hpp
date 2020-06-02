// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Field.hpp
 *
 * \brief This file contains the Abstract base class Field
 *******************************************************************************
 */

#ifndef INLET_FIELD_HPP
#define INLET_FIELD_HPP

#include <string>

namespace axom
{
namespace inlet
{

enum class FieldType {
  Group,
  Int,
  Double,
  String,
  Bool
};

inline std::string fieldTypeToString(FieldType ft)
{
  if (ft == FieldType::Group){
    return "Group";
  }
  else if (ft == FieldType::Int){
    return "Integer";
  }
  else if (ft == FieldType::Double){
    return "Double";
  }
  else if (ft == FieldType::String){
    return "String";
  }
  else if (ft == FieldType::Bool){
    return "Boolean";
  }

  return "Unknown";
}

/*!
 *******************************************************************************
 * \class Field
 *
 * \brief Abstract base class defining the interface of all Map
 *  classes.
 *
 *  Concrete instances need to inherit from this class and implement these
 *  functions.
 *******************************************************************************
 */
class Field
{
public:
  virtual FieldType type() = 0;
  virtual std::string name() = 0;
  virtual std::string description() = 0;
};

} // end namespace inlet
} // end namespace axom

#endif
