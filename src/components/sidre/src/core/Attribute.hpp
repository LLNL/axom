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

private:

  DISABLE_DEFAULT_CTOR(Attribute);
  DISABLE_MOVE_AND_ASSIGNMENT(Attribute);


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
};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_ATTRIBUTE_HPP_ */
