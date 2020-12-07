// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/VariantKey.hpp"

namespace axom
{
namespace inlet
{
VariantKey& VariantKey::operator=(const VariantKey& other)
{
  copyFrom(other);
  return *this;
}
VariantKey& VariantKey::operator=(VariantKey&& other)
{
  copyFrom(std::move(other));
  return *this;
}

VariantKey& VariantKey::operator=(const int key)
{
  if(m_type == VariantKeyType::String)
  {
    using std::string;
    m_string.~string();
  }
  m_int = key;
  return *this;
}

VariantKey& VariantKey::operator=(const std::string& key)
{
  if(m_type == VariantKeyType::String)
  {
    m_string = key;
  }
  else
  {
    new(&m_string) std::string(key);
  }
  return *this;
}

VariantKey& VariantKey::operator=(std::string&& key)
{
  if(m_type == VariantKeyType::String)
  {
    m_string = key;
  }
  else
  {
    new(&m_string) std::string(std::move(key));
  }
  return *this;
}

VariantKey& VariantKey::operator=(const char key[])
{
  if(m_type == VariantKeyType::String)
  {
    m_string = key;
  }
  else
  {
    new(&m_string) std::string(key);
  }
  return *this;
}

VariantKey::~VariantKey()
{
  // Have to call the destructor "manually" as the compiler cannot
  // infer the active member at destruction time
  if(m_type == VariantKeyType::String)
  {
    // Use std::string to get the full basic_string<char, Allocator>
    using std::string;
    m_string.~string();
  }
}

VariantKey::operator int() const
{
  if(m_type != VariantKeyType::Integer)
  {
    SLIC_ERROR(
      "[Inlet] Attempted to retrieve an integer from a non-integer key");
  }
  return m_int;
}

VariantKey::operator const std::string &() const
{
  if(m_type != VariantKeyType::String)
  {
    SLIC_ERROR("[Inlet] Attempted to retrieve a string from a non-string key");
  }
  return m_string;
}

InletType VariantKey::type() const
{
  if(m_type == VariantKeyType::Integer)
  {
    return InletType::Integer;
  }
  else if(m_type == VariantKeyType::String)
  {
    return InletType::String;
  }
  SLIC_ERROR("[Inlet] VariantKey tagged union is in invalid state");
  return InletType::Nothing;
}

bool VariantKey::operator==(const VariantKey& other) const
{
  if(m_type != other.m_type)
  {
    return false;
  }
  if(m_type == VariantKeyType::Integer)
  {
    return m_int == other.m_int;
  }
  else if(m_type == VariantKeyType::String)
  {
    return m_string == other.m_string;
  }
  SLIC_ERROR("[Inlet] VariantKey tagged union is in invalid state");
  return false;
}

void VariantKey::copyFrom(const VariantKey& other)
{
  m_type = other.m_type;
  if(m_type == VariantKeyType::String)
  {
    new(&m_string) std::string(other.m_string);
  }
  else if(m_type == VariantKeyType::Integer)
  {
    m_int = other.m_int;
  }
  else
  {
    SLIC_ERROR("[Inlet] VariantKey tagged union is in invalid state");
  }
}

void VariantKey::copyFrom(VariantKey&& other)
{
  m_type = other.m_type;
  if(m_type == VariantKeyType::String)
  {
    new(&m_string) std::string(std::move(other.m_string));
  }
  else if(m_type == VariantKeyType::Integer)
  {
    m_int = other.m_int;
  }
  else
  {
    SLIC_ERROR("[Inlet] VariantKey tagged union is in invalid state");
  }
}

std::ostream& operator<<(std::ostream& out, const VariantKey& key)
{
  if(key.type() == axom::inlet::InletType::Integer)
  {
    out << static_cast<int>(key);
  }
  else
  {
    out << static_cast<std::string>(key);
  }
  return out;
}

}  // end namespace inlet
}  // end namespace axom
