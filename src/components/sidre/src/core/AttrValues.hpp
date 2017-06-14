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
 * \file AttrValue.hpp
 *
 * \brief   Header file containing definition of AttrValue class.
 *
 ******************************************************************************
 */

// Standard C++ headers
#include <string>
#include <vector>

// Other axom headers
#include "axom/config.hpp"
#include "axom/Macros.hpp"
#include "slic/slic.hpp"

// Sidre project headers
#include "sidre/Attribute.hpp"
#include "sidre/SidreTypes.hpp"

#ifndef SIDRE_ATTRVALUES_HPP_
#define SIDRE_ATTRVALUES_HPP_

namespace axom
{
namespace sidre
{

/*!
 * \class AttrValue
 *
 * \brief Store Attribute values.
 *
 * Each attribute is defined by an instance of Attribute in the DataStore.
 * The attribute has an associated type and index.
 *
 * The assumption is made that attributes will be looked up more often than
 * set.  So getAttribute should be optimized over setAttribute.
 *
 * Another assumption is that the View should not pay for attributes if they
 * are not being used.
 *
 * Space can be minimized by creating more common Attribute first so that they
 * will have a lower index.
 */
class AttrValues
{
public:

  /*!
   * Friend declarations to constrain usage via controlled access to
   * private members.
   */
  friend class View;

private:

  //DISABLE_DEFAULT_CTOR(AttrValues);
  DISABLE_MOVE_AND_ASSIGNMENT(AttrValues);

  /*!
   * \brief Return true if the attribute has been explicitly set; else false.
   */
  bool hasValue(const Attribute * attr) const;

  /*!
   * \brief Create a Node to store attribute.
   */
  bool createNode(const Attribute * attr);

  /*!
   * \brief Set attribute value.
   */
  template<typename ScalarType>
  bool setScalar(const Attribute * attr, ScalarType value)
  {
    DataTypeId arg_id = detail::SidreTT<ScalarType>::id;
    if (arg_id != attr->getDefault().dtype().id())
    {
      SLIC_CHECK_MSG(arg_id == attr->getDefault().dtype().id(),
		     "Incorrect type for attribute '" << attr->getName()
		     << "' of type " << attr->getDefault().dtype().name()
		     << ": " << DataType::id_to_name(arg_id) << ".");
      return false;
    }

    bool ok = createNode(attr);
    if (ok)
    {
      IndexType iattr = attr->getIndex();
      (*m_values)[iattr] = value;
    }
    return ok;
  }

  /*!
   * \brief Set attribute value.
   */
  bool setString(const Attribute * attr, const std::string & value)
  {
    DataTypeId arg_id = CHAR8_STR_ID;
    if (arg_id != attr->getDefault().dtype().id())
    {
      SLIC_CHECK_MSG(arg_id == attr->getDefault().dtype().id(),
		     "Incorrect type for attribute '" << attr->getName()
		     << "' of type " << attr->getDefault().dtype().name()
		     << ": " << DataType::id_to_name(arg_id) << ".");
      return false;
    }

    bool ok = createNode(attr);
    if (ok)
    {
      IndexType iattr = attr->getIndex();
      (*m_values)[iattr] = value;
    }
    return ok;
  }

  /*!
   * \brief Return a value.
   */
  Node::ConstValue getScalar( const Attribute * attr ) const;

  /*!
   * \brief Return a string value.
   */
  const char * getString( const Attribute * attr ) const;

  /*!
   * \brief Return reference to value Node.
   */
  const Node & getValueNodeRef( const Attribute * attr ) const;

//@{
//!  @name Private AttrValues ctor and dtor
//!        (callable only by DataStore methods).

  /*!
   *  \brief Private ctor.
   */
  AttrValues( );

  /*!
   * \brief Private copy ctor.
   */
  AttrValues(const AttrValues& source);

  /*!
   * \brief Private dtor.
   */
  ~AttrValues();

//@}

  ///////////////////////////////////////////////////////////////////
  //
  typedef std::vector< Node > Values;
  ///////////////////////////////////////////////////////////////////

  /// Attributes values.
  Values * m_values;

};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_ATTRVALUES_HPP_ */
