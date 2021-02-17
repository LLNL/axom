// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/Proxy.hpp"

namespace axom
{
namespace inlet
{
InletType Proxy::type() const
{
  // If it's a table, it must be either an object or an array
  if(m_table != nullptr)
  {
    // This is how Inlet stores array types in the datastore
    if(m_table->contains(detail::COLLECTION_GROUP_NAME))
    {
      return InletType::Collection;
    }
    return InletType::Object;
  }
  // Then check if it's a field
  if(m_field != nullptr)
  {
    return m_field->type();
  }
  // Otherwise must be a function
  SLIC_ERROR_IF(!m_func, "[Inlet] Cannot retrieve the type of an empty Proxy");
  return InletType::Function;
}

bool Proxy::contains(const std::string& name) const
{
  if(m_table == nullptr)
  {
    SLIC_ERROR("[Inlet] Cannot index a proxy that refers to a field");
  }
  return m_table->contains(name);
}

Proxy Proxy::operator[](const std::string& name) const
{
  if(m_table == nullptr)
  {
    SLIC_ERROR("[Inlet] Cannot index a proxy that refers to a field");
  }
  return (*m_table)[name];
}

const axom::sidre::Group* Proxy::sidreGroup() const
{
  if(m_table != nullptr)
  {
    return m_table->sidreGroup();
  }
  else if(m_field != nullptr)
  {
    return m_field->sidreGroup();
  }
  else if(m_func != nullptr)
  {
    SLIC_ERROR("[Inlet] Cannot retrieve the sidre::Group for a Function");
  }
  else
  {
    SLIC_ERROR("[Inlet] Cannot retrieve the sidre::Group of an empty Proxy");
  }
  return nullptr;
}

std::string Proxy::name() const
{
  // FIXME: With C++14 we can implement a visit() method that takes a generic lambda
  if(m_table != nullptr)
  {
    return m_table->name();
  }
  else if(m_field != nullptr)
  {
    return m_field->name();
  }
  else if(m_func != nullptr)
  {
    return m_func->name();
  }
  else
  {
    SLIC_ERROR("[Inlet] Cannot retrieve the name of an empty Proxy");
  }
  return "";
}

}  // end namespace inlet
}  // end namespace axom
