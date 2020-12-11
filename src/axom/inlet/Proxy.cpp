// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "Proxy.hpp"

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
    if(m_table->hasTable(detail::CONTAINER_GROUP_NAME))
    {
      return InletType::Container;
    }
    return InletType::Object;
  }
  // Otherwise it must be a field
  if(m_field == nullptr)
  {
    SLIC_ERROR("[Inlet] Cannot retrieve the type of an empty Proxy");
  }
  return m_field->type();
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

}  // end namespace inlet
}  // end namespace axom
