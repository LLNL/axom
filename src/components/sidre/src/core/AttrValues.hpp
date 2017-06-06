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
  bool hasAttrValue(const Attribute * attr) const;

  /*!
   * \brief Set Attribute.
   */
  bool setAttrValue(const Attribute * attr, const std::string & value);

  /*!
   * \brief Return a attribute.
   */
  const char * getAttribute( const Attribute * attr ) const;

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
