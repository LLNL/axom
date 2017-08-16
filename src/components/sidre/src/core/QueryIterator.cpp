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

// Standard C++ headers
//#include <string>
#include <vector>

// Other axom headers
//#include "axom/config.hpp"
//#include "axom/Macros.hpp"
//#include "slic/slic.hpp"

// Sidre project headers
#include "sidre/Group.hpp"
#include "sidre/View.hpp"
#include "sidre/SidreTypes.hpp"

namespace axom
{
namespace sidre
{

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
struct QueryIterator::Cursor
{
  Group *grp;
  IndexType igroup;
  IndexType iview;
  bool grp_done;
  bool first_view;   // true if started processing views
};

struct QueryIterator::CursorStack {
  std::vector< Cursor > stack;
};


/*
 *************************************************************************
 * Find the deepest right-most group.
 *
 * Iterate over the right-most Group until no more Groups are found.
 * Upon return, top of the stack will represent the first View of the
 * deepest Group.  If the Group is empty then iview==InvalidIndex and
 * first_view==true.
 *
 * Any intermediate Groups pushed onto the stack will not have
 * first_view set since their Groups will be visited before their
 * Views.
 *************************************************************************
 */
void QueryIterator::findDeepestGroup(Group *grp)
{
  while (true)
  {
    QueryIterator::Cursor state;
    state.grp    = grp;
    state.igroup = grp->getFirstValidGroupIndex();
    state.iview  = grp->getFirstValidViewIndex();

#if 1
    if (state.iview == InvalidIndex)
    {
      state.grp_done = true;
    }
    else
    {
      state.grp_done = false;
    }
#else
    state.grp_done = false;
#endif
    state.first_view = false;

    m_stack->stack.push_back(state);

    if (state.igroup == InvalidIndex)
    {
      // This Group has no Groups.
      break;
    }

    // Must go deeper...
    grp = grp->getGroup(state.igroup);
  }

  //  if (m_stack->stack.back().iview != InvalidIndex)
  //{
    // The current state represents the first view in a Group.
    m_stack->stack.back().first_view = true;
    //}

}

/*
 *************************************************************************
 * Construct a QueryIterator object.
 *
 * Setup for a depth first query by following down the first group to the bottom.
 *************************************************************************
 */
QueryIterator::QueryIterator(Group * root) :
    m_root(root)
{
  m_stack = new CursorStack;
  findDeepestGroup(root);
}

QueryIterator::~QueryIterator()
{
  delete m_stack;
}


/*
 *************************************************************************
 *  Return true if there is another node to query.
 *************************************************************************
 */
bool QueryIterator::isValid()
{
  if (m_stack->stack.empty())
  {
    return false;
  }
  return true;
}

/*
 *************************************************************************
 * Update iterator to the next available Group or View.
 *
 * If nothing on the stack, must be done.
 * First, look for next Group.  If one is found, start with the deepest node.
 * Once all Groups are exhausted, loop over Views.
 *************************************************************************
 */
void QueryIterator::getNext()
{
  while ( ! m_stack->stack.empty())
  {
    IndexType igroup = m_stack->stack.back().igroup;
    if (igroup != InvalidIndex)
    {
      Group * grp = m_stack->stack.back().grp;
      IndexType inext = grp->getNextValidGroupIndex(igroup);
      m_stack->stack.back().igroup = inext;
      if (inext != InvalidIndex)
      {
	Group *nextgrp = grp->getGroup(inext);
	findDeepestGroup(nextgrp);
	return;    // Found a Group
      }
    }

    IndexType iview = m_stack->stack.back().iview;
    if (iview != InvalidIndex)
    {
      if (m_stack->stack.back().first_view == false) 
      {
	// Finished visiting intermediate groups, now starting on views.
	m_stack->stack.back().first_view = true;
	return;      // iview is the first view
      }
      else
      {
	Group * grp = m_stack->stack.back().grp;
	IndexType inext = grp->getNextValidViewIndex(iview);
	m_stack->stack.back().iview = inext;
	if (inext != InvalidIndex)
        {
	  return;    // Found a View.
	}
      }
    }

    if (m_stack->stack.back().grp_done == false) 
    {
      m_stack->stack.back().grp_done = true;
      return;      // Use Group in this stack entry
    }


#if 0
    // Finished with this Group, go up the tree.
    m_stack->stack.pop_back();
    if (m_stack->stack.empty())
    {
      return;     // No more nodes
    }
#else
    m_stack->stack.pop_back();
#endif


    // Finished with this Group, go up the tree.
    //    m_stack->stack.pop_back();
  }

  return;
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
  if (m_stack->stack.empty())
  {
    return InvalidName;
  }

  QueryIterator::Cursor state = m_stack->stack.back();

  if (state.iview != InvalidIndex)
  {
    View * view = state.grp->getView(state.iview);
    return view->getName();
  }

#if 0
  if (state.igroup != InvalidIndex)
  {
    Group * grp = state.grp->getGroup(state.igroup);
    return grp->getName();
  }
#endif

  return state.grp->getName();
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
