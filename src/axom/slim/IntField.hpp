// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file IntField.hpp
 *
 * \brief This file contains the class IntField
 *******************************************************************************
 */

#ifndef SLIM_INTFIELD_HPP
#define SLIM_INTFIELD_HPP

#include "axom/slim/Field.hpp"

#include <string>

namespace axom
{
namespace slim
{

class IntField
{
public:
  IntField(const std::string& name, const std::string& description)
  : m_name(name)
  , m_description(description) {}

  IntField(std::string&& rname, std::string&& rdescription)
  : m_name(std::move(rname))
  , m_description(std::move(rdescription)) {}

  FieldType type() { return FieldType::Int; }
  std::string name() { return m_name; }
  std::string description() { return m_description; }

  int value() { return m_value; }
  void value(int v) { m_value = v; }

private:
  std::string m_name;
  std::string m_description;

  int m_value;
};

} // end namespace slim
} // end namespace axom

#endif
