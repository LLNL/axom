// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file QueryIterator.hpp
 *
 * \brief   Header file containing definition of QueryIterator class.
 *
 * The QueryIterator is used to visit Views and Groups in a Group.
 *
 * Creating or deleting Views or Groups while iterating over the tree
 * may invalidate the QueryIterator and is not recommended.
 *
 ******************************************************************************
 */

#ifndef SIDRE_ITERATOR_HPP_
#define SIDRE_ITERATOR_HPP_

// Standard C++ headers
#include <string>
#include <stack>

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"

namespace axom
{
namespace sidre
{
class Group;
class View;

/*!
 * \class Iterator
 *
 * \brief Iterate recursively over a Group.
 *
 */
class Iterator
{
public:
  //@{
  //!  @name Iterator ctor and dtor

  /*!
   *  \brief ctor that creates a Iterator which starts at grp.
   */
  Iterator(Group* grp);

  /*!
   * \brief dtor.
   */
  ~Iterator();
  //@}

  //@{
  //!  @name Methods to iterate through a Group.

  /*!
 *  \brief Return true if the Iterator references a Group or View.
 *         Return false if the Iterator has finished its traversal.
 */
  bool isValid();

  /*!
   *  \brief Advance Iterator to the next Group or View.
   *
   *  For each group, its child Groups are traversed first, then its
   *  child Views.  They are traversed in index order.  Typically,
   *  this is the order in which Groups and Views are created.
   *  However, if a Group or View is deleted and a new one created,
   *  the index number will be reused.
   *
   *  After calling advanceToNext, isValid will return true if the
   *  Iterator represents a Group or View and false if the traversal
   *  has finished.  Additional calls to advanceToNext will do nothing
   *  once isValid returns false.
   */
  void advanceToNext();
  //@}

  //@{
  //!  @name Methods to query iterator.

  /*!
   *  \brief Return true if the Iterator references a Group.
   */
  bool isGroup() const;

  /*!
   *  \brief Return true if the Iterator references a View.
   */
  bool isView() const;

  /*!
   *  \brief Return pointer to non-const Group at current iterator position.
   *
   *  If the current position is not a Group, return nullptr.
   */
  Group* asGroup();

  /*!
   *  \brief Return pointer to const Group at current iterator position.
   *
   *  If the current position is not a Group, return nullptr.
   */
  Group const* asGroup() const;

  /*!
   *  \brief Return pointer to non-const View at current iterator position.
   *
   *  If the current position is not a View, return nullptr.
   */
  View* asView();

  /*!
   *  \brief Return pointer to const View at current iterator position.
   *
   *  If the current position is not a View, return nullptr.
   */
  View const* asView() const;

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
  //@}

private:
  DISABLE_DEFAULT_CTOR(Iterator);
  DISABLE_MOVE_AND_ASSIGNMENT(Iterator);

  /*!
   * \brief Push grp on the stack and push child Groups until a Group
   *  is encountered which contains no other Groups.
   */
  void findDeepestGroup(Group* grp);

  ///////////////////////////////////////////////////////////////////
  //
  struct Cursor;
  ///////////////////////////////////////////////////////////////////

  /// Current position in tree
  std::stack<Cursor*> m_stack;
};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_ITERATOR_HPP_ */
