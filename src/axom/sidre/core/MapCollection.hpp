// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file MapCollection.hpp
 *
 * \brief   Header file for MapCollection.
 *
 *          MapCollection is an implemenation of ItemCollection to
 *          hold a collection of items of a fixed type that can be accessed
 *          accessed by string name or sidre::IndexType.
 *
 *          The primary intent is to decouple the implementation of the
 *          collection of times from the Group class which owns collections of
 *          View and child Group objects. This may have other uses,
 *          so it is not dependent on the Group class. This class is
 *          templated on the item type so that it can be used
 *          to hold either View or Group object pointers without
 *          having to code a separate class for each.
 *
 *          \attention This class should be robust against any potential
 *                     user interaction. It doesn't report errors and leaves
 *                     checking of return values to calling code.
 *
 *          \attention Template parameter type must provide a method
 *                     "getName()" that returns a reference to a string object.
 *
 *          \attention The interface is as follows:
 *
 *          \verbatim
 *
 *          - // Return number of items in collection.
 *
 *               size_t getNumItems() const;
 *
 *          - // Return first valid item index (i.e., smallest index over
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
 *          - // Return pointer to item with given name (nullptr if none).
 *
 *               TYPE* getItem(const std::string& name);
 *               TYPE const* getItem(const std::string& name) const ;
 *
 *          - // Return pointer to item with given index (nullptr if none).
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
 *            // (sidre::InvalidIndex if none).
 *
 *               IndexType getItemIndex(const std::string& name) const;
 *
 *          - // Insert item with given name; return index if insertion
 *            // succeeded, and InvalidIndex otherwise.
 *
 *               IndexType insertItem(TYPE* item, const std::string& name);
 *
 *          - // Remove item with given name if it exists and return a
 *            // pointer to it. If it doesn't exist, return nullptr.
 *
 *               TYPE* removeItem(const std::string& name);
 *
 *          - // Remove item with given index if it exists and return a
 *            // pointer to it. If it doesn't exist, return nullptr.
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

#ifndef SIDRE_MAP_COLLECTIONS_HPP_
#define SIDRE_MAP_COLLECTIONS_HPP_

// Standard C++ headers
#include <map>
#include <stack>
#include <string>
#include <vector>

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Types.hpp"

// Sidre project headers
#include "SidreTypes.hpp"
#include "ItemCollection.hpp"

#if defined(AXOM_USE_SPARSEHASH)
  #include <sparsehash/dense_hash_map>
#else
  #include <unordered_map>
#endif

namespace axom
{
namespace sidre
{
////////////////////////////////////////////////////////////////////////
//
// MapCollection keeps an index constant for each item
// as long as it remains in the collection; i.e., don't shift indices
// around.  It has the additional benefit that users can hold on to
// item indices without them being changed without notice.
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
template <typename TYPE>
class MapCollection : public ItemCollection<TYPE>
{
public:
  //
  // Default compiler-generated ctor, dtor, copy ctor, and copy assignment
  // operator suffice for this class.
  //

  ///
  size_t getNumItems() const { return m_items.size() - m_free_ids.size(); }

  ///
  IndexType getFirstValidIndex() const;

  ///
  IndexType getNextValidIndex(IndexType idx) const;

  ///
  bool hasItem(const std::string& name) const
  {
    typename MapType::const_iterator mit = m_name2idx_map.find(name);
    return (mit != m_name2idx_map.end() ? true : false);
  }

  ///
  bool hasItem(IndexType idx) const
  {
    return (idx >= 0 && static_cast<unsigned>(idx) < m_items.size() &&
            m_items[static_cast<unsigned>(idx)]);
  }

  ///
  TYPE* getItem(const std::string& name)
  {
    typename MapType::iterator mit = m_name2idx_map.find(name);
    return (mit != m_name2idx_map.end() ? m_items[mit->second] : nullptr);
  }

  ///
  TYPE const* getItem(const std::string& name) const
  {
    typename MapType::const_iterator mit = m_name2idx_map.find(name);
    return (mit != m_name2idx_map.end() ? m_items[mit->second] : nullptr);
  }

  ///
  TYPE* getItem(IndexType idx)
  {
    return (hasItem(idx) ? m_items[static_cast<unsigned>(idx)] : nullptr);
  }

  ///
  TYPE const* getItem(IndexType idx) const
  {
    return (hasItem(idx) ? m_items[static_cast<unsigned>(idx)] : nullptr);
  }

  ///
  const std::string& getItemName(IndexType idx) const
  {
    return (hasItem(idx) ? m_items[static_cast<unsigned>(idx)]->getName()
                         : InvalidName);
  }

  ///
  IndexType getItemIndex(const std::string& name) const
  {
    typename MapType::const_iterator mit = m_name2idx_map.find(name);
    return (mit != m_name2idx_map.end() ? mit->second : InvalidIndex);
  }

  ///
  IndexType insertItem(TYPE* item, const std::string& name);

  ///
  TYPE* removeItem(const std::string& name);

  ///
  TYPE* removeItem(IndexType idx);

  ///
  void removeAllItems()
  {
    m_items.clear();
    while(!m_free_ids.empty())
    {
      m_free_ids.pop();
    }
#if defined(AXOM_USE_SPARSEHASH)
    if(m_name2idx_map.empty() && m_empty_key != "DENSE_MAP_EMPTY_KEY")
    {
      m_empty_key = "DENSE_MAP_EMPTY_KEY";
      m_name2idx_map.set_empty_key(m_empty_key);
      m_name2idx_map.set_deleted_key("DENSE_MAP_DELETED_KEY");
    }
#endif

    m_name2idx_map.clear();
  }

private:
  std::vector<TYPE*> m_items;
  std::stack<IndexType> m_free_ids;

#if defined(AXOM_USE_SPARSEHASH)
  using MapType = google::dense_hash_map<std::string, IndexType>;
#else
  using MapType = std::unordered_map<std::string, IndexType>;
#endif

  MapType m_name2idx_map;
#if defined(AXOM_USE_SPARSEHASH)
  std::string m_empty_key;
#endif
};

template <typename TYPE>
IndexType MapCollection<TYPE>::getFirstValidIndex() const
{
  IndexType idx = 0;
  while(static_cast<unsigned>(idx) < m_items.size() &&
        m_items[static_cast<unsigned>(idx)] == nullptr)
  {
    idx++;
  }
  return ((static_cast<unsigned>(idx) < m_items.size()) ? idx : InvalidIndex);
}

template <typename TYPE>
IndexType MapCollection<TYPE>::getNextValidIndex(IndexType idx) const
{
  if(idx == InvalidIndex)
  {
    return InvalidIndex;
  }

  idx++;
  while(static_cast<unsigned>(idx) < m_items.size() &&
        m_items[static_cast<unsigned>(idx)] == nullptr)
  {
    idx++;
  }
  return ((static_cast<unsigned>(idx) < m_items.size()) ? idx : InvalidIndex);
}

template <typename TYPE>
IndexType MapCollection<TYPE>::insertItem(TYPE* item, const std::string& name)
{
  bool use_recycled_index = false;
  IndexType idx = m_items.size();
  if(!m_free_ids.empty())
  {
    idx = m_free_ids.top();
    m_free_ids.pop();
    use_recycled_index = true;
  }

#if defined(AXOM_USE_SPARSEHASH)
  if(m_name2idx_map.empty() && m_empty_key != "DENSE_MAP_EMPTY_KEY")
  {
    m_empty_key = "DENSE_MAP_EMPTY_KEY";
    m_name2idx_map.set_empty_key(m_empty_key);
    m_name2idx_map.set_deleted_key("DENSE_MAP_DELETED_KEY");
  }
#endif

  if(m_name2idx_map.insert(std::make_pair(name, idx)).second)
  {
    // name was inserted into map
    if(use_recycled_index)
    {
      m_items[idx] = item;
    }
    else
    {
      m_items.push_back(item);
    }
    return idx;
  }
  else
  {
    // name was NOT inserted into map, return free index if necessary
    if(use_recycled_index)
    {
      m_free_ids.push(idx);
    }
    return InvalidIndex;
  }
}

template <typename TYPE>
TYPE* MapCollection<TYPE>::removeItem(const std::string& name)
{
  TYPE* ret_val = nullptr;

  typename MapType::iterator mit = m_name2idx_map.find(name);
  if(mit != m_name2idx_map.end())
  {
    IndexType idx = mit->second;

    ret_val = m_items[idx];

    m_name2idx_map.erase(mit);
    m_items[idx] = nullptr;
    m_free_ids.push(idx);
  }

  return ret_val;
}

template <typename TYPE>
TYPE* MapCollection<TYPE>::removeItem(IndexType idx)
{
  if(hasItem(idx))
  {
    TYPE* item = removeItem(m_items[idx]->getName());
    return item;
  }
  else
  {
    return nullptr;
  }
}

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_MAP_COLLECTIONS_HPP_ */
