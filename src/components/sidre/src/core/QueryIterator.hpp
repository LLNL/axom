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
 * \file QueryIterator.hpp
 *
 * \brief   Header file containing definition of QueryIterator class.
 *
 * The QueryIterator is used to visit Views and Groups in a Group.
 *
 ******************************************************************************
 */

#ifndef SIDRE_QUERYITERATOR_HPP_
#define SIDRE_QUERYITERATOR_HPP_

// Standard C++ headers
//#include <string>

// Other axom headers
#include "axom/config.hpp"
#include "axom/Macros.hpp"
//#include "slic/slic.hpp"

// Sidre project headers
//#include "sidre/SidreTypes.hpp"

namespace axom
{
namespace sidre
{

class Group;

/*!
 * \class Attribute
 *
 * \brief Attributes for a View
 *
 */
class QueryIterator
{
public:

private:
  DISABLE_DEFAULT_CTOR(QueryIterator);
  DISABLE_MOVE_AND_ASSIGNMENT(QueryIterator);

  /// start of iteration
  Group * root;

//@{
//!  @name Private QueryIterator ctor and dtor
//!        (callable only by Group methods).

  /*!
   *  \brief Private ctor that creates an QueryIterator which starts at grp.
   */
  QueryIterator( Group * grp);

  /*!
   * \brief Private dtor.
   */
  ~QueryIterator();
//@}

};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_QUERYITERATOR_HPP_ */
