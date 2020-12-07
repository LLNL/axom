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

#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \enum InletType
 *
 * \brief Enumeration of basic types for things in inlet
 *******************************************************************************
 */
enum class InletType
{
  Nothing,
  Bool,
  String,
  Integer,
  // TODO: Unsigned integer
  Double,
  Object,
  Container
};

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
  /*!
   *****************************************************************************
   * \brief Parameterized constructors for initializing the variant
   * \param [in] key The key to initialize with 
   *****************************************************************************
   */
  VariantKey(const int key) : m_int(key), m_type(VariantKeyType::Integer) { }
  /// \overload
  VariantKey(const std::string& key)
    : m_string(key)
    , m_type(VariantKeyType::String)
  { }
  /// \overload
  VariantKey(std::string&& key)
    : m_string(std::move(key))
    , m_type(VariantKeyType::String)
  { }
  /// \overload
  // Explicit const char[] overload required to convert to std::string instead of ptr -> int
  VariantKey(const char key[]) : m_string(key), m_type(VariantKeyType::String)
  { }

  VariantKey(const VariantKey& other) { copyFrom(other); }
  VariantKey(VariantKey&& other) { copyFrom(std::move(other)); }

  /*!
   *****************************************************************************
   * \brief Copy assignment operator
   *****************************************************************************
   */
  VariantKey& operator=(const VariantKey& other);
  /*!
   *****************************************************************************
   * \brief Move assignment operator
   *****************************************************************************
   */
  VariantKey& operator=(VariantKey&& other);

  VariantKey& operator=(const int key);

  VariantKey& operator=(const std::string& key);
  VariantKey& operator=(std::string&& key);

  VariantKey& operator=(const char key[]);
  ~VariantKey();

  operator int() const;
  operator const std::string &() const;

  InletType type() const;

  bool operator==(const VariantKey& other) const;

private:
  /*!
   *****************************************************************************
   * \brief Move assignment operator
   *****************************************************************************
   */
  void copyFrom(const VariantKey& other);
  void copyFrom(VariantKey&& other);

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

std::ostream& operator<<(std::ostream& out, const VariantKey& key);

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
