// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Associated header file
#include "Iterator.hpp"

// Sidre project headers
#include "axom/sidre/core/Group.hpp"
#include "axom/sidre/core/View.hpp"
#include "axom/sidre/core/SidreTypes.hpp"

namespace axom
{
namespace sidre
{
/*
 * Groups are traversed in a depth first order by keeping a stack of
 * Groups in the current path starting with the initial group at the
 * bottom of the stack and the deepest Group at the top.
 *
 * m_igroup is the current Group being visited. When a Group is
 * visited, a new Cursor is pushed onto the stack with m_grp pointing
 * to the Group.  After the top of the stack is fully traversed, it is
 * popped from the stack then the new top of stack will have its
 * m_igroup incremented via getNextValidGroupIndex.  Each group will
 * eventually be at the top of the stack.
 *
 * Views are traversed starting at the initial View found by
 * getFirstValidViewIndex then traversed using getNextValidViewIndex.
 * View visits do not push anything additional onto the Cursor stack.
 *
 * m_iview is the current View being visited.  If m_iview is
 * InvalidIndex, then the Cursor represents the Group m_grp.
 *
 * m_is_first_view will be set true to indicate that the first View
 * has been visited.  This will be true for the deepest Group returned
 * from findDeepestGroup if it contains Views.  Any intermediate
 * Groups between the start Group and the deepest Group will have
 * m_is_first_view set to false. This is checked when an intermediate
 * Group becomes the top of the stack, its Groups are first traversed
 * then its Views.  m_is_first_view is set to true by advanceToNext to
 * flag that the first View has been visited and the next call to
 * advanceToNext should update iview before using it.
 *
 * Once all of the Views for a group are visited, m_is_group is set to
 * true to indicate that Cursor represents the parent Group.  The next
 * call to advanceToNext will then pop the stack and update igroup.
 * m_is_group is also set to true by findDeepestGroup if the Group
 * found has no Views to signal advanceToNext that the Group has been
 * visited.
 */
struct Iterator::Cursor
{
  Group* m_grp;
  IndexType m_igroup;
  IndexType m_iview;
  bool m_is_group;
  bool m_is_first_view;

  /*
 *************************************************************************
 * Return true if Cursor references a Group.
 *************************************************************************
 */
  bool isGroup() { return m_iview == InvalidIndex ? true : false; }

  /*
 *************************************************************************
 * Return true if Cursor references a View.
 *************************************************************************
 */
  bool isView() { return m_iview == InvalidIndex ? false : true; }

  /*
 *************************************************************************
 *  If the Cursor represents a Group, return the Group. Else,
 *  return nullptr
 *************************************************************************
 */
  Group* getCurrentGroup() { return m_iview == InvalidIndex ? m_grp : nullptr; }

  /*
 *************************************************************************
 *  If the Cursor represents a View, return the View. Else,
 *  return nullptr
 *************************************************************************
 */
  View* getCurrentView()
  {
    return m_iview == InvalidIndex ? nullptr : m_grp->getView(m_iview);
  }
};

/*
 *************************************************************************
 * Find the next Group in a depth-first traversal of grp.
 *
 * If the Group has Views, then m_is_first_view is set to true to indicate
 * that the Cursor represents the first view in a Group. If it does
 * not have Views, then m_is_group is set to true indicating that the
 * Cursor represents the Group.  i.e. it has no Groups or Views.
 *
 *************************************************************************
 */
void Iterator::findDeepestGroup(Group* grp)
{
  while(true)
  {
    Iterator::Cursor* state = new Iterator::Cursor;
    state->m_grp = grp;
    state->m_igroup = grp->getFirstValidGroupIndex();
    state->m_iview = grp->getFirstValidViewIndex();
    state->m_is_group = false;
    state->m_is_first_view = false;

    m_stack.push(state);

    if(state->m_igroup == InvalidIndex)
    {
      // This Group has no Groups.
      break;
    }

    // Must go deeper...
    grp = grp->getGroup(state->m_igroup);
  }

  // Check last Group pushed to see if it represents a Group or View
  Iterator::Cursor* state = m_stack.top();
  if(state->isView())
  {
    // Start by visiting the first View
    state->m_is_first_view = true;
  }
  else
  {
    state->m_is_group = true;
  }
}

/*
 *************************************************************************
 * Construct a Iterator instance.
 *
 * Setup for a depth first query by following down grp to the bottom.
 *************************************************************************
 */
Iterator::Iterator(Group* grp) { findDeepestGroup(grp); }

/*
 *************************************************************************
 * Destructor for Iterator.
 *************************************************************************
 */
Iterator::~Iterator()
{
  // If the Iterator is destoryed before it has iterated over all
  // Groups or Views, m_stack will not be empty.
  while(!m_stack.empty())
  {
    delete m_stack.top();
    m_stack.pop();
  }
}

/*
 *************************************************************************
 *  Return true if the Iterator references a Group or View.
 *  Return false if the Iterator has finished its traversal.
 *************************************************************************
 */
bool Iterator::isValid() { return !m_stack.empty(); }

/*
 *************************************************************************
 * Advance iterator to the next Group or View.
 *
 * If nothing on the stack, nothing to do.  Else, look for next Group.
 * If one is found, start with its deepest Group or View.  Once all
 * Groups are exhausted, loop over Views. Finally, the Cursor
 * represents the Group since all of its child Groups and Views have
 * been visited.
 *************************************************************************
 */
void Iterator::advanceToNext()
{
  while(!m_stack.empty())
  {
    Iterator::Cursor* state = m_stack.top();
    Group* grp = state->m_grp;

    if(state->m_igroup != InvalidIndex)
    {
      state->m_igroup = grp->getNextValidGroupIndex(state->m_igroup);
      if(state->m_igroup != InvalidIndex)
      {
        Group* nextgrp = grp->getGroup(state->m_igroup);
        findDeepestGroup(nextgrp);
        return;  // Found a Group
      }
    }

    if(state->m_iview != InvalidIndex)
    {
      if(state->m_is_first_view == false)
      {
        state->m_is_first_view = true;
        return;  // iview is the first view
      }
      else
      {
        state->m_iview = grp->getNextValidViewIndex(state->m_iview);
        if(state->m_iview != InvalidIndex)
        {
          return;  // Found a View.
        }
      }
    }

    if(state->m_is_group == false)
    {
      state->m_is_group = true;
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
 * Return true if Iterator references a Group.
 *
 *************************************************************************
 */
bool Iterator::isGroup() const
{
  if(m_stack.empty())
  {
    return false;
  }

  return m_stack.top()->isGroup();
}

/*
 *************************************************************************
 *
 * Return true if Iterator references a View.
 *
 *************************************************************************
 */
bool Iterator::isView() const
{
  if(m_stack.empty())
  {
    return false;
  }

  return m_stack.top()->isView();
}

/*
 *************************************************************************
 *
 *  Return pointer to non-const Group at current iterator position.
 *
 *************************************************************************
 */
Group* Iterator::asGroup()
{
  if(m_stack.empty())
  {
    return nullptr;
  }

  return m_stack.top()->getCurrentGroup();
}

/*
 *************************************************************************
 *
 *  Return pointer to const Group at current iterator position.
 *
 *************************************************************************
 */
Group const* Iterator::asGroup() const
{
  if(m_stack.empty())
  {
    return nullptr;
  }

  return m_stack.top()->getCurrentGroup();
}

/*
 *************************************************************************
 *
 *  Return pointer to non-const View at current iterator position.
 *
 *************************************************************************
 */
View* Iterator::asView()
{
  if(m_stack.empty())
  {
    return nullptr;
  }

  return m_stack.top()->getCurrentView();
}

/*
 *************************************************************************
 *
 *  Return pointer to const View at current iterator position.
 *
 *************************************************************************
 */
View const* Iterator::asView() const
{
  if(m_stack.empty())
  {
    return nullptr;
  }

  return m_stack.top()->getCurrentView();
}

/*
 *************************************************************************
 *
 * Return const reference to name of current iterator object.
 *
 *************************************************************************
 */
const std::string& Iterator::getName() const
{
  if(m_stack.empty())
  {
    return InvalidName;
  }

  Iterator::Cursor* state = m_stack.top();

  if(state->isView())
  {
    View* view = state->m_grp->getView(state->m_iview);
    return view->getName();
  }

  return state->m_grp->getName();
}

#if 0
/*
 *************************************************************************
 *
 * Return path of Group object, not including its name.
 *
 *************************************************************************
 */
std::string Iterator::getPath() const
{
  #if 1
  std::string thePath("here");
  #else
  const Group* root = getDataStore()->getRoot();
  const Group* curr = getParent();
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
