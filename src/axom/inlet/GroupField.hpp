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

#ifndef INLET_GROUPFIELD_HPP
#define INLET_GROUPFIELD_HPP

#include "axom/inlet/Field.hpp"

#include <string>

namespace axom
{
namespace inlet
{

class GroupField : public Field
{
public:
  GroupField(const std::string& name, const std::string& description)
    : m_name(name)
    , m_description(description) {}

  GroupField(std::string&& rname, std::string&& rdescription)
    : m_name(std::move(rname))
    , m_description(std::move(rdescription)) {}

  FieldType type() { return FieldType::Int; }
  std::string name() { return m_name; }
  std::string description() { return m_description; }

private:
  std::string m_name;
  std::string m_description;
};

} // end namespace inlet
} // end namespace axom

#endif
