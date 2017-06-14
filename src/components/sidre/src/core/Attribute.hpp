/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 ******************************************************************************
 *
 * \file Attribute.hpp
 *
 * \brief   Header file containing definition of Attribute class.
 *
 ******************************************************************************
 */

// Standard C++ headers
#include <string>

// Other axom headers
#include "axom/config.hpp"
#include "axom/Macros.hpp"
#include "slic/slic.hpp"

// Sidre project headers
#include "sidre/SidreTypes.hpp"

#ifndef SIDRE_ATTRIBUTE_HPP_
#define SIDRE_ATTRIBUTE_HPP_

namespace axom
{
namespace sidre
{

/*!
 * \class Attribute
 *
 * \brief Attributes for a View
 *
 */
class Attribute
{
public:

  /*!
   * Friend declarations to constrain usage via controlled access to
   * private members.
   */
  friend class DataStore;

  /*!
   * \brief Return const reference to name of Attribute object.
   */
  const std::string& getName() const
  {
    return m_name;
  }

  /*!
   * \brief Return the unique index of this Attribute object.
   */
  IndexType getIndex() const
  {
    return m_index;
  }

  /*!
   * \brief Set default value of Attribute. Return true if successfully changed.
   *
   * The type of the default cannot be changed after the attribute is created.
   */
  template<typename ScalarType>
  bool setDefault(ScalarType value)
  {
    DataTypeId arg_id = detail::SidreTT<ScalarType>::id;
    if (m_default_value.dtype().is_empty() || 
        arg_id == m_default_value.dtype().id())
    {
      m_default_value = value;
      return true;
    }
    else
    {
      SLIC_CHECK_MSG(arg_id == m_default_value.dtype().id(),
		     "Cannot change type of attribute '" << m_name
		     << "' from " << m_default_value.dtype().name()
		     << " to " << DataType::id_to_name(arg_id) << ".");
      return false;
    }
  }

  /*!
   * \brief Set default value of Attribute. Return true if successfully changed.
   *
   * The type of the default cannot be changed after the attribute is created.
   */
  bool setDefault(const std::string & value)
  {
    DataTypeId arg_id = CHAR8_STR_ID;
    if (m_default_value.dtype().is_empty() || 
        arg_id == m_default_value.dtype().id())
    {
      m_default_value = value;
      return true;
    }
    else
    {
      SLIC_CHECK_MSG(arg_id == m_default_value.dtype().id(),
		     "Cannot change type of attribute '" << m_name
		     << "' from " << m_default_value.dtype().name()
		     << " to " << DataType::id_to_name(arg_id) << ".");
      return false;
    }
  }

  /*!
   * \brief Return default value of Attribute.
   */
  const Node & getDefault() const
  {
    return m_default_value;
  }

private:

  DISABLE_DEFAULT_CTOR(Attribute);
  DISABLE_MOVE_AND_ASSIGNMENT(Attribute);

  /*!
   * \brief Set index of attribute within DataStore.
   *        Called as part of DataStore->createAttribute.
   */
  void setIndex(IndexType index)
  {
    m_index = index;
  }

//@{
//!  @name Private Attribute ctor and dtor
//!        (callable only by DataStore methods).

  /*!
   *  \brief Private ctor that creates an Attribute with given name
   *         which has no data associated with it.
   */
  Attribute( const std::string& name );

  /*!
   * \brief Private copy ctor.
   */
  Attribute(const Attribute& source);

  /*!
   * \brief Private dtor.
   */
  ~Attribute();

//@}

  /// Name of this Attribute object.
  std::string m_name;

  /// Attribute's unique index within DataStore object that created it.
  IndexType m_index;

  /// Default value of attribute.
  Node m_default_value;

};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_ATTRIBUTE_HPP_ */
