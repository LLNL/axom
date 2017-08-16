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
 * Creating Views or Groups while traversing the tree may invalidate
 * the QueryIterator and is not recommended.
 *
 ******************************************************************************
 */

#ifndef SIDRE_QUERYITERATOR_HPP_
#define SIDRE_QUERYITERATOR_HPP_

// Standard C++ headers
#include <string>
#include <vector>

// Other axom headers
#include "axom/config.hpp"
#include "axom/Macros.hpp"
//#include "slic/slic.hpp"

// Sidre project headers
#include "sidre/SidreTypes.hpp"

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

//@{
//!  @name Private QueryIterator ctor and dtor
//!        (callable only by Group methods).

  /*!
   *  \brief Private ctor that creates an QueryIterator which starts at grp.
   */
  QueryIterator(Group * grp);

#if 0
  /*!
   * \brief Private dtor.
   */
  ~QueryIterator();
#endif
//@}

  /*!
   *  \brief Return true if the Query references a Group or View.
   *         Return false if there are no more nodes to visit.
   */
  bool isValid();

  /*!
   *  \brief Update iterator to the next available node.
   */
  void getNext();


  /*!
   * \brief Return const reference to name of current iterator object.
   *
   * \sa getPath(), getPathName()
   */
  const std::string& getName() const;

#if 0
  /*!
   * \brief Return path of query object, not including its name.
   *
   * \sa getName(), getPathName()
   */
  std::string getPath() const;

  /*!
   * \brief Return full path of Group object, including its name.
   *
   * If a DataStore contains a Group tree structure a/b/c/d/e, the 
   * following results are expected:
   *
   * Method Call      | Result
   * -----------------|----------
   * e->getName()     | e
   * e->getPath()     | a/b/c/d
   * e->getPathName() | a/b/c/d/e
   *
   * \sa getName(), getPath(), View::getPathName()
   */
  std::string getPathName() const
  {
    if (getPath().length() < 1)
    {
      return getName();
    }

    return getPath() + getPathDelimiter() + getName();
  }
#endif

private:
  DISABLE_DEFAULT_CTOR(QueryIterator);
  DISABLE_MOVE_AND_ASSIGNMENT(QueryIterator);

  /*!
   * \brief Push grp on the stack and find the deepest right-most group.
   */
  void findDeepestGroup(Group *grp);

  /*
   *  first_view will be true if the Curor represents the first view
   *  of a Group.  This will be true for the deepest Group returned
   *  from findDeepestGroup.  Note that iview may be InvalidGroup if
   *  the group has no views.  Intermediate Groups will process their
   *  own Groups first before proceeding to processing their Views.
   *
   *  If false, the first View has already been visited and it is
   *  necessary to call getNextValidViewIndex.
   *
   *  finished_views will be true when there are no views to process.
   *  Either the group initially had no views, or all of the views
   *  have been visited.  The Cursor represent the Group grp, and not
   *  any views in the Group.
   */
  struct Cursor
  {
    Group *grp;
    IndexType igroup;
    IndexType iview;
    bool grp_done;
    bool first_view;   // true if started processing views
  };

  /// start of iteration
  Group * m_root;

  ///////////////////////////////////////////////////////////////////
  //
  typedef std::vector< Cursor > StateStack;
  ///////////////////////////////////////////////////////////////////

  /// Current position in tree
  StateStack m_stack; 

};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_QUERYITERATOR_HPP_ */
