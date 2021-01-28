// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Proxy.hpp
 *
 * \brief This file contains the class definition of Inlet's Proxy class.
 *******************************************************************************
 */

#ifndef INLET_PROXY_HPP
#define INLET_PROXY_HPP

#include <type_traits>

#include "axom/inlet/Field.hpp"
#include "axom/inlet/Table.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class Proxy
 *
 * \brief Provides a uniform interface for access and conversion to primitive
 * and user-defined types
 *
 * \see Inlet Field
 * \see Inlet Table
 *******************************************************************************
 */
class Proxy
{
public:
  Proxy() = default;
  /*!
   *******************************************************************************
   * \brief Constructs a proxy view onto a table
   * 
   * \param [in] table The table to construct a proxy into
   *******************************************************************************
   */
  Proxy(Table& table) : m_table(&table) { }

  /*!
   *******************************************************************************
   * \brief Constructs a proxy view onto a field
   * 
   * \param [in] field The field to construct a proxy into
   *******************************************************************************
   */
  Proxy(Field& field) : m_field(&field) { }

  /*!
   *******************************************************************************
   * \brief Constructs a proxy view onto a function
   * 
   * \param [in] func The function to construct a proxy into
   *******************************************************************************
   */
  Proxy(Function& func) : m_func(&func) { }

  /*!
   *******************************************************************************
   * \brief Returns an object from the proxy
   * 
   * \tparam T The type of the object to retrieve
   * \return The retrieved object
   * \note Implicit conversions are defined only for primitive types - to 
   * convert to a user-defined type, use an explicit conversion with static_cast
   * or T{}
   *******************************************************************************
   */
  template <
    typename T,
    typename SFINAE = typename std::enable_if<
      detail::is_inlet_primitive<T>::value || detail::is_inlet_primitive_array<T>::value ||
      detail::is_inlet_primitive_dict<T>::value || detail::is_std_function<T>::value ||
      detail::is_primitive_std_vector<T>::value>::type>
  operator T() const
  {
    return get<T>();
  }

  /*!
   *****************************************************************************
   * \brief Return whether a subobject with the given name is present in 
   * the table referred to by the calling proxy.
   *
   * \return Boolean value indicating whether this Table's subtree contains a
   * Field or Table with the given name.
   *****************************************************************************
   */
  bool contains(const std::string& name) const;

  /*!
   *****************************************************************************
   * \brief Returns the type of the stored value
   * 
   * \return The type
   * \see InletType
   *****************************************************************************
   */
  InletType type() const;

  /*!
   *******************************************************************************
   * \brief Obtains a proxy view into the proxy for either a Field/Table subobject
   * 
   * Returns a reference via a lightweight proxy object to the element in the 
   * datastore at the index specified by the name.  This can be a field 
   * or a table.
   * 
   * \param [in] name The name of the subobject
   * \return A view onto the subobject
   *******************************************************************************
   */
  Proxy operator[](const std::string& name) const;

  /*!
   *******************************************************************************
   * \brief Returns a user-defined type from the proxy
   * 
   * \tparam T The type of the object to retrieve
   * \return The retrieved object
   * \pre The Proxy must refer to a table object
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<!detail::is_inlet_primitive<T>::value &&
                            !detail::is_std_function<T>::value,
                          T>::type
  get() const
  {
    SLIC_ASSERT_MSG(m_table != nullptr,
                    "[Inlet] Tried to read a user-defined type from a Proxy "
                    "containing a single field or function");
    return m_table->get<T>();
  }

  /*!
   *******************************************************************************
   * \brief Returns a function type from the proxy
   * 
   * \tparam T The type of the object to retrieve
   * \return The retrieved object
   * \pre The Proxy must refer to a function object
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<detail::is_std_function<T>::value, T>::type get() const
  {
    SLIC_ASSERT_MSG(m_func != nullptr,
                    "[Inlet] Tried to read a function from a Proxy "
                    "containing a field or table");
    return m_func->get<typename detail::std_function_signature<T>::type>();
  }

  /*!
   *******************************************************************************
   * \brief Calls the function
   * 
   * \param [in] args The parameter pack for the function's arguments
   * \tparam Ret The user-specified return type, needed to fully disambiguate the
   * function to call
   * \tparam Args The types of the user-specified arguments, deduced automatically
   * 
   * \return The function's result
   *******************************************************************************
   */
  template <typename Ret, typename... Args>
  Ret call(Args&&... args) const
  {
    SLIC_ASSERT_MSG(m_func != nullptr,
                    "[Inlet] Tried to call a Proxy "
                    "containing a field or table");
    return m_func->call<Ret>(std::forward<Args>(args)...);
  }

  /*!
   *******************************************************************************
   * \brief Returns a primitive type from the proxy
   * 
   * \tparam T The type of the object to retrieve
   * \return The retrieved object
   * \pre The Proxy must refer to a field object
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<detail::is_inlet_primitive<T>::value, T>::type get() const
  {
    SLIC_ASSERT_MSG(m_field != nullptr,
                    "[Inlet] Tried to read a primitive type from a Proxy "
                    "containing a table or function");
    return m_field->get<T>();
  }

private:
  Table* m_table = nullptr;
  Field* m_field = nullptr;
  Function* m_func = nullptr;
};

}  // end namespace inlet
}  // end namespace axom

#endif  // INLET_PROXY_HPP
