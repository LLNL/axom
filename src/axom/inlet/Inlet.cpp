// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Inlet.cpp
 *
 * \brief This file contains the class implementation of Inlet, the main class
 *        for the Inlet component.
 *******************************************************************************
 */

#include "axom/inlet/Inlet.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{

GroupField* Inlet::addGroup(const std::string& name,
                            const std::string& description)
{
  SLIC_ASSERT_MSG(m_map != nullptr, "Map not set");
  SLIC_ASSERT_MSG(m_backend != nullptr, "Backend not set");

  GroupField* group = new GroupField(name, description);
  m_backend->add((Field*)group);
  return group;
}

GroupField* Inlet::addGroup(std::string&& rname,
                            std::string&& rdescription)
{
  SLIC_ASSERT_MSG(m_map != nullptr, "Map not set");
  SLIC_ASSERT_MSG(m_backend != nullptr, "Backend not set");

  GroupField* group = new GroupField(rname, rdescription);
  m_backend->add((Field*)group);
  return group;
}

IntField* Inlet::addIntField(const std::string& name,
                             const std::string& description,
                             int defaultValue)
{
  SLIC_ASSERT_MSG(m_map != nullptr, "Map not set");
  SLIC_ASSERT_MSG(m_backend != nullptr, "Backend not set");

  IntField* intField = new IntField(name, description);
  int v;
  if(m_map->getInt(name, v))
  {
    intField->value(v);
  }
  else
  {
    intField->value(defaultValue);
  }
  m_backend->add((Field*)intField);
  return intField;
}

IntField* Inlet::addIntField(const std::string& name,
                             const std::string& description,
                             bool required)
{
  SLIC_ASSERT_MSG(m_map != nullptr, "Map not set");
  SLIC_ASSERT_MSG(m_backend != nullptr, "Backend not set");

  IntField* intField = nullptr;
  int v;
  if(m_map->getInt(name, v))
  {
    intField = new IntField(name, description);
    intField->value(v);
    m_backend->add((Field*)intField);
  }
  else if(required)
  {
    SLIC_ERROR("Required field is no found in input deck: " + name);
  }
  return intField;
}

IntField* Inlet::getIntField(const std::string& name)
{
  SLIC_ASSERT_MSG(m_backend != nullptr, "Backend not set");

  Field* field = m_backend->get(name);
  if (field == nullptr)
  {
    return nullptr;
  }
  if (field->type() != FieldType::Int)
  {
    std::string error_msg = fmt::format("Integer Field named '{0}' was asked for"
                                        " but recieved type {1}",
                                        name,
                                        fieldTypeToString(field->type()));
    SLIC_ERROR(error_msg);
  }

  return (IntField*)field;
}

} // end namespace inlet
} // end namespace axom
