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
 * Each attribute is defined by an instance if Attribute in the DataStore.
 * The attribute has an associated type and index.
 *
 * While attributes are associated with a View, they are saved in the Group.
 * Each attribute has a std::vector of values indexed by the View's index.
 * Attributes are saved in a std::vector of vectors so that attributes.
 */
class AttrValues
{
public:

  /*!
   * Friend declarations to constrain usage via controlled access to
   * private members.
   */
  friend class Group;

private:

  //DISABLE_DEFAULT_CTOR(AttrValues);
  DISABLE_MOVE_AND_ASSIGNMENT(AttrValues);

  bool setAttrValue(const Attribute * attr, IndexType idx, const std::string & value);

  const std::string & getAttribute( const Attribute * attr, IndexType idx ) const;

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
  typedef std::vector< std::vector<std::string> * > Attributes;
  ///////////////////////////////////////////////////////////////////

  /// Name of this Attribute object.
  Attributes * m_attributes;

};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_ATTRVALUES_HPP_ */
