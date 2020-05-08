// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Structure.cpp
 *
 * \brief This file contains the class implementation of Structure.
 *******************************************************************************
 */

#include "axom/slim/Structure.hpp"

#include "axom/slic.hpp"

namespace axom
{
namespace slim
{

GroupField* Structure::addGroup(const std::string& name,
                                            const std::string& description)
{
  GroupField* group = new GroupField(name, description);
  m_fields.push_back((Field*)group);
  return group;
}

GroupField* Structure::addGroup(std::string&& rname,
                                            std::string&& rdescription)
{
  GroupField* group = new GroupField(rname, rdescription);
  m_fields.push_back((Field*)group);
  return group;
}

IntField* Structure::addIntField(const std::string& name,
                                 const std::string& description,
                                 int defaultValue)
{
  SLIC_ASSERT_MSG(m_map != nullptr, "Map not set");

  IntField* intField = new IntField(name, description);
  int v;
  if(m_map->getInt(name, v)) {
    intField->value(v);  
  }
  else {
    intField->value(defaultValue);
  }
  m_fields.push_back((Field*)intField);
  return intField;
}

IntField* Structure::addIntField(const std::string& name,
                                 const std::string& description)
{
  SLIC_ASSERT_MSG(m_map != nullptr, "Map not set");

  IntField* intField = nullptr;
  int v;
  if(m_map->getInt(name, v)) {
    intField = new IntField(name, description);
    intField->value(v);  
    m_fields.push_back((Field*)intField);
  }
  return intField;
}

} // end namespace slim
} // end namespace axom
