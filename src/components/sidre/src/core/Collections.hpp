/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
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
 * \file
 *
 * \brief   Header file for Collection classes.
 *
 *          Each of these classes holds a collection of items of a fixed
 *          type that can be accessed by string name or sidre::IndexType.
 *
 *          The primary intent is to decouple the implementation of the
 *          collections from the DataGroup class which owns collections of
 *          DataView and child DataGroup objects. They may have other uses,
 *          so they are not dependent on the DataGroup class. Each class is
 *          templated on the item type so that the same class can be used
 *          to hold either DataView or DataGroup object pointers without
 *          having to code a separate class for each.
 *
 *          By having various collections that obey the same interface,
 *          we can explore alternative collection implementations for
 *          performance (insertion, lookup, etc.) and memory overhead.
 *          The collection used by the DataGroup class can be changed via
 *          the collection typedef in the DataGroup class header file.
 *
 *          To try another collection, encapsulate it in a new class with
 *          the API described below or pass it as a template parameter to
 *          an existing class below if that works.
 *
 *          IMPORTANT: These classes should be robust against any potential
 *                     user interaction. They don't report errors and leave
 *                     checking of return values to calling code.
 *
 *          IMPORTANT: Template parameter type must provide a method
 *                     "getName()" that returns a reference to a string object.
 *
 *          IMPORTANT: The common interface each collection class provides
 *                     is as follows:
 *
 *          \verbatim
 *
 *          - // Return number of items in collection.
 *
 *               size_t getNumItems() const;
 *
 *          - // Retturn first valid item index (i.e., smallest index over
 *            // all items).
 *            // sidre::InvalidIndex returned if no items in collection
 *
 *               IndexType getFirstValidIndex() const;
 *
 *          - // Return next valid item index after given index (i.e., smallest
 *            // index over all indices larger than given one).
 *            // sidre::InvalidIndex returned
 *
 *               IndexType getNextValidIndex(IndexType idx) const;
 *
 *          - // Return true if item with given name in collection; else false.
 *
 *               bool hasItem(const std::string& name) const;
 *
 *          - // Return true if item with given index in collection; else false.
 *
 *               bool hasItem(IndexType idx) const;
 *
 *          - // Return pointer to item with given name (ATK_NULLPTR if none).
 *
 *               TYPE* getItem(const std::string& name);
 *               TYPE const* getItem(const std::string& name) const ;
 *
 *          - // Return pointer to item with given index (ATK_NULLPTR if none).
 *
 *               TYPE* getItem(IndexType idx);
 *               TYPE const* getItem(IndexType idx) const;
 *
 *          - // Return name of object with given index
 *            // (sidre::InvalidName if none).
 *
 *               std::string getItemName(IndexType idx) const;
 *
 *          - // Return index of object with given name
 *            // (sidre::InvalidName if none).
 *
 *               IndexType getItemIndex(const std::string& name) const;
 *
 *          - // Insert item with given name; return true if insertion
 *            // succeeded, and false otherwise.
 *
 *               bool insertItem(TYPE* item, const std::string& name);
 *
 *          - // Remove item with given name if it exists and return a
 *            // pointer to it. If it doesn't exist, return ATK_NULLPTR.
 *
 *               TYPE* removeItem(const std::string& name);
 *
 *          - // Remove item with given name if it exists and return a
 *            // pointer to it. If it doesn't exist, return ATK_NULLPTR.
 *
 *               TYPE* removeItem(IndexType idx);
 *
 *          - // Remove all items (items not destroyed).
 *
 *               void removeAllItems();
 *
 *          - // Clear all items and destroy them.
 *
 *               void deleteAllItems();
 *
 *          \endverbatim
 *
 ******************************************************************************
 */

#ifndef COLLECTIONS_HPP_
#define COLLECTIONS_HPP_

// Standard C++ headers
#include <map>
#include <stack>
#include <string>
#include <vector>

// other CS Toolkit headers
#include "common/CommonTypes.hpp"

// SiDRe project headers
#include "SidreTypes.hpp"


namespace asctoolkit
{
namespace sidre
{

////////////////////////////////////////////////////////////////////////
//
// MapCollection implementation follows original implementation
// in DataGroup.
//
////////////////////////////////////////////////////////////////////////

/*!
 *************************************************************************
 *
 * \class MapCollection
 *
 * \brief MapCollection is a container class template for holding
 *        a collection of items of template parameter type TYPE, using
 *        a map container of type MAP_TYPE.
 *
 * \warning Only std::map and std::unordered_map have been tried so far.
 *          These classes have identical APIs for the functionality we
 *          are using.
 *
 *************************************************************************
 */
template <typename TYPE, typename MAP_TYPE>
class MapCollection
{
public:

  //
  // Default compiler-generated ctor, dtor, copy ctor, and copy assignment
  // operator suffice for this class.
  //

  ///
  size_t getNumItems() const
  {
    return m_items.size();
  }

  ///
  IndexType getFirstValidIndex() const
  {
    return getNextValidIndex(-1);
  }

  ///
  IndexType getNextValidIndex(IndexType idx) const
  {
    return hasItem(idx+1) ? idx+1 : InvalidIndex;
  }

  ///
  bool hasItem(const std::string& name) const
  {
    typename MAP_TYPE::const_iterator mit = m_name2idx_map.find(name);
    return ( mit != m_name2idx_map.end() ? true : false );
  }

  ///
  bool hasItem(IndexType idx) const
  {
    return (idx >= 0 &&
            static_cast<unsigned>(idx) < getNumItems() &&
            m_items[idx]);
  }

  ///
  TYPE * getItem(const std::string& name)
  {
    typename MAP_TYPE::iterator mit = m_name2idx_map.find(name);
    return ( mit != m_name2idx_map.end() ?
             m_items[ mit->second ] : ATK_NULLPTR );
  }

  ///
  TYPE const * getItem(const std::string& name) const
  {
    typename MAP_TYPE::const_iterator mit = m_name2idx_map.find(name);
    return ( mit != m_name2idx_map.end() ?
             m_items[ mit->second ] : ATK_NULLPTR );
  }

  ///
  TYPE * getItem(IndexType idx)
  {
    return ( hasItem(idx) ? m_items[idx] : ATK_NULLPTR );
  }

  ///
  TYPE const * getItem(IndexType idx) const
  {
    return ( hasItem(idx) ? m_items[idx] : ATK_NULLPTR );
  }

  ///
  const std::string& getItemName(IndexType idx) const
  {
    return ( hasItem(idx) ? m_items[idx]->getName() : InvalidName );
  }

  ///
  IndexType getItemIndex(const std::string& name) const
  {
    typename MAP_TYPE::const_iterator mit = m_name2idx_map.find(name);
    return ( mit != m_name2idx_map.end() ?
             mit->second : InvalidIndex );
  }

  ///
  bool insertItem(TYPE * item, const std::string& name);

  ///
  TYPE * removeItem(const std::string& name);

  ///
  TYPE * removeItem(IndexType idx);

  ///
  void removeAllItems()
  {
    m_items.clear();
    m_name2idx_map.clear();
  }

private:
  std::vector<TYPE *>  m_items;
  MAP_TYPE m_name2idx_map;
};

template <typename TYPE, typename MAP_TYPE>
bool MapCollection<TYPE, MAP_TYPE>::insertItem(TYPE * item,
                                               const std::string& name)
{
  if ( m_name2idx_map.insert( std::make_pair(name, m_items.size()) ).second )
  {
    // name was inserted into map
    m_items.push_back( item );
    return true;
  }
  else
  {
    // name was NOT inserted into map
    return false;
  }
}

template <typename TYPE, typename MAP_TYPE>
TYPE * MapCollection<TYPE, MAP_TYPE>::removeItem(const std::string& name)
{
  TYPE * ret_val = ATK_NULLPTR;

  typename MAP_TYPE::iterator mit = m_name2idx_map.find(name);
  if ( mit != m_name2idx_map.end() )
  {
    IndexType idx = mit->second;

    ret_val = m_items[idx];

    m_name2idx_map.erase(mit);
    m_items.erase(m_items.begin() + idx);

    // Decrement approriate item indices
    for (mit = m_name2idx_map.begin() ; mit != m_name2idx_map.end() ; ++mit)
    {
      if (mit->second > idx)
      {
        mit->second--;
      }
    }
  }

  return ret_val;
}

template <typename TYPE, typename MAP_TYPE>
TYPE * MapCollection<TYPE, MAP_TYPE>::removeItem(IndexType idx)
{
  if ( hasItem(idx) )
  {
    TYPE * item = removeItem( m_items[idx]->getName() );
    return item;
  }
  else
  {
    return ATK_NULLPTR;
  }
}



////////////////////////////////////////////////////////////////////////
//
// NewMapCollection is an attempt to improve performance over
// MapCollection above by keeping an index constant for each item
// as long as it iremains in the collection; i.e., don't shift indices
// around.  It has the additional benefit that users can hold on to
// item indices without them being changed without notice.
//
////////////////////////////////////////////////////////////////////////

/*!
 *************************************************************************
 *
 * \class NewMapCollection
 *
 * \brief NewMapCollection is a container class template for holding
 *        a collection of items of template parameter type TYPE, using
 *        a map container of type MAP_TYPE.
 *
 * \warning Only std::map and std::unordered_map have been tried so far.
 *          These classes have identical APIs for the functionality we
 *          are using.
 *
 *************************************************************************
 */
template <typename TYPE, typename MAP_TYPE>
class NewMapCollection
{
public:

  //
  // Default compiler-generated ctor, dtor, copy ctor, and copy assignment
  // operator suffice for this class.
  //

  ///
  size_t getNumItems() const
  {
    return m_items.size() - m_free_ids.size();
  }

  ///
  IndexType getFirstValidIndex() const;

  ///
  IndexType getNextValidIndex(IndexType idx) const;

  ///
  bool hasItem(const std::string& name) const
  {
    typename MAP_TYPE::const_iterator mit = m_name2idx_map.find(name);
    return ( mit != m_name2idx_map.end() ? true : false );
  }

  ///
  bool hasItem(IndexType idx) const
  {
    return (idx >= 0 &&
            static_cast<unsigned>(idx) < m_items.size() &&
            m_items[idx]);
  }

  ///
  TYPE * getItem(const std::string& name)
  {
    typename MAP_TYPE::iterator mit = m_name2idx_map.find(name);
    return ( mit != m_name2idx_map.end() ?
             m_items[ mit->second ] : ATK_NULLPTR );
  }

  ///
  TYPE const * getItem(const std::string& name) const
  {
    typename MAP_TYPE::const_iterator mit = m_name2idx_map.find(name);
    return ( mit != m_name2idx_map.end() ?
             m_items[ mit->second ] : ATK_NULLPTR );
  }

  ///
  TYPE * getItem(IndexType idx)
  {
    return ( hasItem(idx) ? m_items[idx] : ATK_NULLPTR );
  }

  ///
  TYPE const * getItem(IndexType idx) const
  {
    return ( hasItem(idx) ? m_items[idx] : ATK_NULLPTR );
  }

  ///
  const std::string& getItemName(IndexType idx) const
  {
    return ( hasItem(idx) ? m_items[idx]->getName() : InvalidName );
  }

  ///
  IndexType getItemIndex(const std::string& name) const
  {
    typename MAP_TYPE::const_iterator mit = m_name2idx_map.find(name);
    return ( mit != m_name2idx_map.end() ?
             mit->second : InvalidIndex );
  }

  ///
  bool insertItem(TYPE * item, const std::string& name);

  ///
  TYPE * removeItem(const std::string& name);

  ///
  TYPE * removeItem(IndexType idx);

  ///
  void removeAllItems()
  {
    m_items.clear();
    while ( !m_free_ids.empty() )
    {
      m_free_ids.pop();
    }
    m_name2idx_map.clear();
  }

private:
  std::vector<TYPE *>  m_items;
  std::stack< IndexType > m_free_ids;
  MAP_TYPE m_name2idx_map;
};

template <typename TYPE, typename MAP_TYPE>
IndexType NewMapCollection<TYPE, MAP_TYPE>::getFirstValidIndex() const
{
  return getNextValidIndex(-1);
}

template <typename TYPE, typename MAP_TYPE>
IndexType NewMapCollection<TYPE, MAP_TYPE>::getNextValidIndex(IndexType idx) const
{
  idx++;
  while ( static_cast<unsigned>(idx) < m_items.size() &&
          m_items[idx] == ATK_NULLPTR )
  {
    idx++;
  }
  return ( (static_cast<unsigned>(idx) < m_items.size()) ? idx : InvalidIndex );
}


template <typename TYPE, typename MAP_TYPE>
bool NewMapCollection<TYPE, MAP_TYPE>::insertItem(TYPE * item,
                                                  const std::string& name)
{
  bool use_recycled_index = false;
  IndexType idx = m_items.size();
  if ( !m_free_ids.empty() )
  {
    idx = m_free_ids.top();
    m_free_ids.pop();
    use_recycled_index = true;
  }

  if ( m_name2idx_map.insert( std::make_pair(name, idx) ).second )
  {
    // name was inserted into map
    if ( use_recycled_index )
    {
      m_items[idx] = item;
    }
    else
    {
      m_items.push_back(item);
    }
    return true;
  }
  else
  {
    // name was NOT inserted into map, return free index if necessary
    if ( use_recycled_index )
    {
      m_free_ids.push(idx);
    }
    return false;
  }
}

template <typename TYPE, typename MAP_TYPE>
TYPE * NewMapCollection<TYPE, MAP_TYPE>::removeItem(const std::string& name)
{
  TYPE * ret_val = ATK_NULLPTR;

  typename MAP_TYPE::iterator mit = m_name2idx_map.find(name);
  if ( mit != m_name2idx_map.end() )
  {
    IndexType idx = mit->second;

    ret_val = m_items[idx];

    m_name2idx_map.erase(mit);
    m_items[idx] = ATK_NULLPTR;
    m_free_ids.push(idx);
  }

  return ret_val;
}

template <typename TYPE, typename MAP_TYPE>
TYPE * NewMapCollection<TYPE, MAP_TYPE>::removeItem(IndexType idx)
{
  if ( hasItem(idx) )
  {
    TYPE * item = removeItem( m_items[idx]->getName() );
    return item;
  }
  else
  {
    return ATK_NULLPTR;
  }
}

} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* COLLECTIONS_HPP_ */
