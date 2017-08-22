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
 * \file QueryIterator.cpp
 *
 * \brief   Implementation file for QueryIterator class.
 *
 ******************************************************************************
 */

// Associated header file
#include "QueryIterator.hpp"

// Sidre project headers
#include "sidre/Group.hpp"
#include "sidre/View.hpp"
#include "sidre/SidreTypes.hpp"

namespace axom
{
namespace sidre
{

/*
 * Groups are traversed in a depth first order by keeping a stack of
 * Groups in the current path starting with the initial group at the
 * bottom of the stack and the right-most deepest Group at the top.
 * Each group will eventually be at the top of the stack.
 *
 * Views are traversed starting at the initial View found by 
 * getFirstValidViewIndex then traversed using getNextValidViewIndex.
 * View visits do not push anything onto the Cursor stack.
 *
 * iview is the current View being visited.
 * igroup is the current Group being visited.
 *
 * first_view will be set true to indicate that the first View has
 * been visited.  This will be true for the deepest Group returned
 * from findDeepestGroup if it contains Views.  Any intermediate
 * Groups between the initial Group and the deepest Group will have
 * first_view set to false. This is necessary because when any
 * intermediate Group is at the top of the stack, its Group are first
 * traversed, then its Views.  first_view is set to true by getNext to
 * flag that the first View has been visited and the next call to
 * getNext should update iview before using it.
 *
 * Once all of the Views for a group are visited, is_group is set to
 * true to indicate that Cursor represents the parent Group.  The next
 * call to getNext will then pop the stack and update igroup.
 * is_group is also set to true by findDeepestGroup if the Group found
 * has no Views to signal getNext that the Group has been visited.
 */
struct QueryIterator::Cursor
{
  Group *grp;
  IndexType igroup;
  IndexType iview;
  bool is_group;
  bool first_view;
};

/*
 *************************************************************************
 * Find the deepest right-most group.
 *
 * Iterate over the right-most Group until no more Groups are found.
 *
 * If the Group has Views, then first_view is set to true to indicate
 * that the Cursor represents the first view in a Group. If it does
 * not have Views, then is_group is set to true indicating that the
 * Cursor represents the Group.  i.e. it has no Groups or Views.
 *
 *************************************************************************
 */
void QueryIterator::findDeepestGroup(Group *grp)
{
  while (true)
  {
    QueryIterator::Cursor * state = new QueryIterator::Cursor;
    state->grp    = grp;
    state->igroup = grp->getFirstValidGroupIndex();
    state->iview  = grp->getFirstValidViewIndex();
    state->is_group = false;
    state->first_view = false;

    m_stack.push(state);

    if (state->igroup == InvalidIndex)
    {
      // This Group has no Groups.
      break;
    }

    // Must go deeper...
    grp = grp->getGroup(state->igroup);
  }

  // Check last Group pushed to see if it represents a Group or View
  QueryIterator::Cursor * state = m_stack.top();
  if (state->iview == InvalidIndex)
  {
    // No Views in the group
    state->is_group = true;
  }
  else
  {
    // Start by visiting the first View
    state->first_view = true;
  }

}

/*
 *************************************************************************
 * Construct a QueryIterator instance.
 *
 * Setup for a depth first query by following down root to the bottom.
 *************************************************************************
 */
QueryIterator::QueryIterator(Group * root)
{
  findDeepestGroup(root);
}

/*
 *************************************************************************
 * Destructor for QueryIterator.
 *************************************************************************
 */
QueryIterator::~QueryIterator()
{
  // If the Iterator is destory before it has iterated over all nodes,
  // m_stack will not be empty.
  while ( ! m_stack.empty())
  {
    delete m_stack.top();
    m_stack.pop();
  }
}


/*
 *************************************************************************
 *  Return true if there is another node to query.
 *************************************************************************
 */
bool QueryIterator::isValid()
{
  if (m_stack.empty())
  {
    return false;
  }
  return true;
}

/*
 *************************************************************************
 * Update iterator to the next available Group or View.
 *
 * If nothing on the stack, nothing to do.  Else, look for next Group.
 * If one is found, start with its deepest node.  Once all Groups are
 * exhausted, loop over Views. Finally, the Cursor represents the Group
 * since all of its child Groups and Views have been visited.
 *************************************************************************
 */
void QueryIterator::getNext()
{
  while ( ! m_stack.empty())
  {
    QueryIterator::Cursor * state = m_stack.top();

    Group * grp = state->grp;
    IndexType igroup = state->igroup;
    if (igroup != InvalidIndex)
    {
      IndexType inext = grp->getNextValidGroupIndex(igroup);
      state->igroup = inext;
      if (inext != InvalidIndex)
      {
	Group *nextgrp = grp->getGroup(inext);
	findDeepestGroup(nextgrp);
	return;    // Found a Group
      }
    }

    IndexType iview = state->iview;
    if (iview != InvalidIndex)
    {
      if (state->first_view == false) 
      {
	state->first_view = true;
	return;      // iview is the first view
      }
      else
      {
	IndexType inext = grp->getNextValidViewIndex(iview);
	state->iview = inext;
	if (inext != InvalidIndex)
        {
	  return;    // Found a View.
	}
      }
    }

    if (state->is_group == false) 
    {
      state->is_group = true;
      return;
    }

    delete state;
    m_stack.pop();
  }

  return;
}

/*
 *************************************************************************
 *
 * Return true if QueryIterator references a Group.
 *
 *************************************************************************
 */
bool QueryIterator::isGroup() const
{
  if (m_stack.empty())
  {
    return false;
  }

  if (m_stack.top()->iview != InvalidIndex)
  {
    return false;
  }

  return true;
}

/*
 *************************************************************************
 *
 * Return true if QueryIterator references a View.
 *
 *************************************************************************
 */
bool QueryIterator::isView() const
{
  if (m_stack.empty())
  {
    return false;
  }

  if (m_stack.top()->iview != InvalidIndex)
  {
    return true;
  }

  return false;
}

/*
 *************************************************************************
 *
 *  Return pointer to non-const Group at current iterator position.
 *
 *************************************************************************
 */
Group * QueryIterator::asGroup()
{
  if (m_stack.empty())
  {
    return AXOM_NULLPTR;
  }

  QueryIterator::Cursor * state = m_stack.top();

  if (state->iview != InvalidIndex)
  {
    return AXOM_NULLPTR;
  }

  return state->grp;
}

/*
 *************************************************************************
 *
 *  Return pointer to const Group at current iterator position.
 *
 *************************************************************************
 */
Group const * QueryIterator::asGroup() const
{
  if (m_stack.empty())
  {
    return AXOM_NULLPTR;
  }

  QueryIterator::Cursor * state = m_stack.top();

  if (state->iview != InvalidIndex)
  {
    return AXOM_NULLPTR;
  }

  return state->grp;
}

/*
 *************************************************************************
 *
 *  Return pointer to non-const View at current iterator position.
 *
 *************************************************************************
 */
View * QueryIterator::asView()
{
  if (m_stack.empty())
  {
    return AXOM_NULLPTR;
  }

  QueryIterator::Cursor * state = m_stack.top();

  if (state->iview != InvalidIndex)
  {
    return state->grp->getView(state->iview);
  }

  return AXOM_NULLPTR;
}

/*
 *************************************************************************
 *
 *  Return pointer to const View at current iterator position.
 *
 *************************************************************************
 */
View const * QueryIterator::asView() const
{
  if (m_stack.empty())
  {
    return AXOM_NULLPTR;
  }

  QueryIterator::Cursor * state = m_stack.top();

  if (state->iview != InvalidIndex)
  {
    return state->grp->getView(state->iview);
  }

  return AXOM_NULLPTR;
}

/*
 *************************************************************************
 *
 * Return const reference to name of current iterator object.
 *
 *************************************************************************
 */
const std::string & QueryIterator::getName() const
{
  if (m_stack.empty())
  {
    return InvalidName;
  }

  QueryIterator::Cursor * state = m_stack.top();

  if (state->iview != InvalidIndex)
  {
    View * view = state->grp->getView(state->iview);
    return view->getName();
  }

  return state->grp->getName();
}

#if 0
/*
 *************************************************************************
 *
 * Return path of Group object, not including its name.
 *
 *************************************************************************
 */
std::string QueryIterator::getPath() const
{
#if 1
  std::string thePath("here");
#else
  const Group * root = getDataStore()->getRoot();
  const Group * curr = getParent();
  std::string thePath = curr->getName();
  curr = curr->getParent();

  while (curr != root)
  {
    thePath = curr->getName() + s_path_delimiter + thePath;
    curr = curr->getParent();
  }
#endif

  return thePath;
}
#endif

} /* end namespace sidre */
} /* end namespace axom */
