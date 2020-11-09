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
    if(m_table->hasTable("_inlet_array"))
    {
      return InletType::Array;
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

}  // end namespace inlet
}  // end namespace axom
