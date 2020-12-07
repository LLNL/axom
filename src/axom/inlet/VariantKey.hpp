// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file VariantKey.hpp
 *
 * \brief This file contains the class definition of VariantKey, Inlet's generic
 * associative array key type
 *******************************************************************************
 */

#ifndef INLET_KEY_HPP
#define INLET_KEY_HPP

#include "axom/inlet/Field.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class VariantKey
 *
 * \brief Provides a generic key type for mixed-key associative arrays/dictionaries
 * 
 * \note Can be thought of as limited implementation of std::variant<int, std::string>
 *******************************************************************************
 */
class VariantKey
{
public:
  VariantKey(const int key) : m_int(key), m_type(VariantKeyType::Integer) { }
  VariantKey(const std::string& key)
    : m_string(key)
    , m_type(VariantKeyType::String)
  { }
  VariantKey(std::string&& key)
    : m_string(std::move(key))
    , m_type(VariantKeyType::String)
  { }
  VariantKey(const char key[]) : m_string(key), m_type(VariantKeyType::String)
  { }
  VariantKey(const VariantKey& other) { copy_from(other); }
  VariantKey(VariantKey&& other) { copy_from(std::move(other)); }
  VariantKey& operator=(const VariantKey& other)
  {
    copy_from(other);
    return *this;
  }
  VariantKey& operator=(VariantKey&& other)
  {
    copy_from(std::move(other));
    return *this;
  }

  VariantKey& operator=(const int key)
  {
    if(m_type == VariantKeyType::String)
    {
      using std::string;
      m_string.~string();
    }
    m_int = key;
    return *this;
  }

  VariantKey& operator=(const std::string& key)
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

  VariantKey& operator=(std::string&& key)
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

  VariantKey& operator=(const char key[])
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

  ~VariantKey()
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

  operator int() const
  {
    if(m_type != VariantKeyType::Integer)
    {
      SLIC_ERROR(
        "[Inlet] Attempted to retrieve an integer from a non-integer key");
    }
    return m_int;
  }

  operator const std::string &() const
  {
    if(m_type != VariantKeyType::String)
    {
      SLIC_ERROR(
        "[Inlet] Attempted to retrieve a string from a non-string key");
    }
    return m_string;
  }

  InletType type() const
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

  bool operator==(const VariantKey& other) const
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

private:
  void copy_from(const VariantKey& other)
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

  void copy_from(VariantKey&& other)
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

  // Subset of InletType
  enum class VariantKeyType
  {
    Integer,
    String
  };
  /*!
   *****************************************************************************
   * \brief Member of sum (variant) type for each possible key type
   *****************************************************************************
   */
  union
  {
    int m_int;
    std::string m_string;
  };

  // Active member of the union
  VariantKeyType m_type;
};

}  // end namespace inlet
}  // end namespace axom

namespace std
{
template <>
struct hash<axom::inlet::VariantKey>
{
  size_t operator()(const axom::inlet::VariantKey& key) const
  {
    if(key.type() == axom::inlet::InletType::Integer)
    {
      return std::hash<int> {}(key);
    }
    else
    {
      return std::hash<std::string> {}(key);
    }
  }
};
}  // end namespace std

#endif  // INLET_KEY_HPP
