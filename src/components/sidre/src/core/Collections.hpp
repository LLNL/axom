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
 *          - // Return name of object with given index (InvalidName if none).
 *
 *               std::string getItemName(IndexType idx) const;
 *
 *          - // Return index of object with given name (InvalidIndex if none).
 *
 *               IndexType getItemIndex(const std::string& name) const;
 *
 *          - // Insert item with given name; return true if insertion
 *               succeeded, and false otherwise.
 *
 *               bool insertItem(TYPE* item, const std::string& name);
 *
 *          - // Remove item with given name if it exists and return a
 *               pointer to it. If it doesn't exist, return ATK_NULLPTR.
 *
 *               TYPE* removeItem(const std::string& name);
 *
 *          - // Remove item with given name if it exists and return a
 *               pointer to it. If it doesn't exist, return ATK_NULLPTR.
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
// MapCollection class definition and implementation.
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
    // item was inserted into map
    m_items.push_back( item );
    return true;
  }
  else
  {
    // item was NOT inserted into map
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
// Other implementations forthcoming....
//
////////////////////////////////////////////////////////////////////////


} /* end namespace sidre */
} /* end namespace asctoolkit */

#endif /* COLLECTIONS_HPP_ */
