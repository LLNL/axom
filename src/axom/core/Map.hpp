// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MAP_HPP_
#define AXOM_MAP_HPP_

#include "axom/config.hpp"                  // for compile-time defines
#include "axom/core/Macros.hpp"             // for axom macros
#include "axom/core/memory_management.hpp"  // for memory allocation functions
#include "axom/core/Types.hpp"              // for axom types

// C/C++ includes
#include <functional>  //for hashing until a custom method is defined
#include <iostream>

namespace axom
{
namespace experimental
{
namespace axom_map
{
/// \name Map Supporting Data Structures
/// @{
/*!
   * \class Node
   * 
   * \brief Node is the foundational data type of the Map, being the lowest-level data storage unit,
   *  with the key and value data for a given item, along with a next value for the linked list.
   *
   * \tparam Key the type of key for the key-value pair stored within the Node instance.
   * \tparam T the type of value for the key-value pair stored within the Node instance.
   */
template <typename Key, typename T>
struct Node
{
  IndexType next;
  Key key;
  T value;
  bool operator==(const Node& rhs)
  {
    return (this->next == rhs.next) && (this->key == rhs.key) &&
      (this->value == rhs.value);
  }
};

/*!
   * \class Pair
   * 
   * \brief Pair is strictly for returning from insert, to match the pair-returning setup of STL 
   *  unordered_map.
   * 
   * \tparam Key the type of key stored in the associated Node.
   * \tparam T the type of value stored in the associated Node.
   */
template <typename Key, typename T>
struct Pair
{
  axom_map::Node<Key, T>* first;
  bool second;

  Pair(axom_map::Node<Key, T>* node, bool status) : first(node), second(status)
  { }
};

/*!
   * \class Bucket
   *
   * \brief Bucket is the implementation of a bucket for this Map with chaining. It's a singly-linked
   * list, and supports search, insertion, and removal.
   *
   * \tparam Key the type of the keys stored in each Node in the list.
   * \tparapm T the type of the values stored in each Node in the list.
   */
template <typename Key, typename T>
class Bucket
{
public:
  /// \name Bucket Construction
  ///@{
  /*!
     * \brief Default constructor so new Buckets just construct into nothing. 
     *  Can likely just be removed. 
     */
  Bucket() { }

  /*!
     * \brief Basic constructor, creates singly-linked list with provided size. 
     *
     * \param [in] len the number of items this Bucket instance needs to be able to store.
     *
     * \pre len > 0
     */
  Bucket(int len)
  {
    init(len);
  }

  /*!
     * \brief If a Bucket instance is created with the default contructor, this
     *  is used to initialize and allocate all the member variables. 
     *
     * \param [in] len the number of items this Bucket instance needs to be able to store.
     *
     * \pre len > 0
     */
  void init(int len)
  {
    m_list = axom::allocate<axom_map::Node<Key, T>>(len);
    for(int i = 0; i < len - 1; i++)
    {
      m_list[i].next = i + 1;
    }
    m_list[len - 1].next = -1;
    m_free = 0;
    m_head = -1;
    m_capacity = len;
    m_size = 0;
    m_end.key = Key {0};
    m_end.value = T {0};
    m_end.next = -2;
  }
  /// @}
  /// \name Bucket Methods
  /// @{
  /*!
     * \brief Inserts a given key-value pair into the linked list if an item with 
     *  given key does not already exist. Fails if list is full or item with given key
     *  already exists.
     *
     * \param [in] key the key of the key-value pair to be inserted.
     * \param [in] value the value of the key-value pair to be inserted.
     *
     * \return Returns a Pair object, with the first element being a pointer to the inserted node
     *  and the second element being boolean value true if successful. If insertion fails,
     *  the bool is set to False, and the node pointed to is the end node if the list was full,
     *  or the item occupying the queried slot if an item with the given key already existed.
     */
  axom_map::Pair<Key, T> insert_no_update(Key key, T value)
  {
    if(m_free != -1)
    {
      IndexType ind = m_head;
      if(m_head == -1)
      {
        m_head = m_free;
        ind = m_free;
        m_free = m_list[m_free].next;
        m_list[m_head].next = -1;
        m_list[m_head].key = key;
        m_list[m_head].value = value;
        m_size++;
      }
      else
      {
        while(m_list[ind].next != -1)
        {
          if(m_list[ind].key == key)
          {
            axom_map::Pair<Key, T> ret(&(m_list[ind]), false);
            return ret;
          }
          ind = m_list[ind].next;
        }
        m_list[ind].next = m_free;
        ind = m_free;
        m_free = m_list[m_free].next;
        m_list[ind].next = -1;
        m_list[ind].key = key;
        m_list[ind].value = value;
        m_size++;
      }
      axom_map::Pair<Key, T> ret(&(m_list[ind]), true);
      return ret;
    }
    axom_map::Pair<Key, T> ret(&(m_end), false);
    return ret;
  }

  /*!
     * \brief Removes item with the given key from the linked list, if item exists. Relinks accordingly.
     *
     * \param [in] key the key of the item to be found and removed.
     *
     * \return Returns true if item was found and removed, false otherwise. 
     *  
     */
  bool remove(Key key)
  {
    if(m_head == -1)
    {
      return false;
    }
    IndexType ind = m_head;
    IndexType prev = -1;
    do
    {
      if(m_list[ind].key == key)
      {
        if(prev != -1)
        {
          m_list[prev].next = m_list[ind].next;
        }
        if(ind == m_head)
        {
          m_head = m_list[ind].next;
        }
        m_list[ind].next = m_free;
        m_free = ind;
        m_size--;
        return true;
      }
      prev = ind;
      ind = m_list[ind].next;
    } while(ind != -1);
    return false;
  }

  /*!
    * \brief Returns pointer to node containing key-value pair matching provided key, if possible.
    *
    * \param [in] key the key of the item to be searched for.
    *
    * \return Returns pointer to the requested node if possible, and to special end node otherwise.
    *
    */
  axom_map::Node<Key, T>& find(Key key)
  {
    if(m_head == -1)
    {
      return m_end;
    }
    IndexType ind = m_head;
    do
    {
      if(m_list[ind].key == key)
      {
        return m_list[ind];
      }

      ind = m_list[ind].next;
    } while(ind != -1);
    return m_end;
  }

  /*!
   * \brief Returns maximum number of items that can be stored in this Bucket instance.
   *
   * \return Returns maximum number of items that can be stored in this Bucket instance.
   *
   */
  int get_capacity() { return m_capacity; }

  /*!
   * \brief Returns current number of items in this Bucket instance.
   *
   * \return Returns current number of items in this Bucket instance.
   */
  int get_size() { return m_size; }
  ///@}

  axom_map::Node<Key, T>* m_list;
  axom_map::Node<Key, T> m_end;
  IndexType m_head;
  IndexType m_free;
  int m_capacity, m_size;
};

/// @}
}  // namespace axom_map

/*!
 * \class Map
 *
 * \brief Provides a hashmap implementation, relying upon chaining to 
 *  resolve hash collisions. Simple hashmap for now, future work will use 
 *  RAJA and Umpire to allow use of map within kernels. Resembles 
 *  unordered_map, which we can't use due to a lack of host-device decoration.
 *
 * \tparam Key the type of keys.
 * \tparam T the type of values to hold.
 * \tparam Hash functor that takes an object of Key type and returns 
 *  a hashed value of size_t (for now)
 */

template <typename Key, typename T, typename Hash = std::hash<Key>>
class Map
{
public:
  /// \name Map Constructors
  /// @{

  /*!
   *\brief Constructs a Map instance with the given number of buckets.
   *
   * \param [in] num_buckets the number of buckets used in the Map.
   * \param [in] bucket_len the maximum number of items that can co-exist in a bucket.
   *
   * \pre num_buckets > 0
   * \pre bucket_len > 0
   */
  Map(int num_buckets, int bucket_len = 10)
  {
    init(num_buckets, bucket_len);
  }
  /// @}
  
  /*!
   * Destructor. Frees associated buffers.
   */ 
  ~Map()
  {
    deallocate();
  }

  /// \name Map Methods
  ///@{
  /*!
   * \brief Resizes Map by creating new Map instance and inserting all values
   * from original Map.
   *
   * \param [in] buckets user-specified number of buckets in new Map. -1 to fall back to multiplicative factor.
   * \param [in] factor user-specified multiplicative factor to determine size of new Map based upon existing. 
   *
   * \note if neither input is specified, the new Map will be of size 2*old Map size. 
   * \note both map instances exist in memory until rehash returns. Risks running out of memory.
   */
  void rehash(int buckets = -1, int factor = -1)
  {
    int newlen = 0;
    IndexType ind;
    std::size_t hashed;
    if(buckets != -1 && buckets > m_bucket_count)
    {
      newlen = buckets;
    }
    else if(factor != -1 && factor > 1)
    {
      newlen = factor * m_bucket_count;
    }
    else
    {
      newlen = 2 * m_bucket_count;
    }

    axom_map::Bucket<Key, T>* new_list = alloc_map(newlen, m_bucket_len);

    for(int i = 0; i < m_bucket_count; i++)
    {
      ind = m_buckets[i].m_head;
      while(ind != -1)
      {
        hashed = (get_hash(m_buckets[i].m_list[ind].key) % newlen);
        new_list[hashed].insert_no_update(m_buckets[i].m_list[ind].key,
                                          m_buckets[i].m_list[ind].value);
        ind = m_buckets[i].m_list[ind].next;
      }
      axom::deallocate(m_buckets[i].m_list);
    }
    axom::deallocate(m_buckets);
    m_buckets = new_list;
    m_bucket_count = newlen;
    return;
  }
  /*!
   * \brief Deallocates memory resources for this Map instance, resets state as if constructed with 0 elements.
   *  Also used in destructor.
   *
   */ 
  void deallocate(){
    for(int i = 0; i < m_bucket_count; i++){
      axom::deallocate(m_buckets[i].m_list);
    }
    axom::deallocate(m_buckets);
    m_bucket_count = 0;
    m_bucket_len = 0;
    m_size = 0;
    m_load_factor = 0;
  }

  /*!
   * \brief Initializes this Map instance, mostly in case deallocate was used and there's a desire to continue use.
   *
   */ 
  void init(int num_buckets, int bucket_len = 10){
    m_bucket_count = num_buckets;
    m_bucket_len = bucket_len;
    m_size = 0;
    m_load_factor = 0;
    m_buckets = alloc_map(m_bucket_count, m_bucket_len);
  }
  /*!
   * \brief Inserts given key-value pair into the Map instance.
   *
   * \param [in] key the key used to index into the Map to store this item
   * \param [in] val the data value of the pair to insert into the map
   *
   * \note Deviation from STL unordered_map is that insertion fails when map full, rather than automatically rehashing. 
   *  Manual call of rehash() required.
   *
   * \return A pair whose first element is a pointer to the inserted item and bool set to true
   *  if successful. Otherwise, the second element is set to false, and the first to one of two values: a sentinel node 
   *  if the map is overfilled, and the actual node with the given key if an item with given key already exists.
   *
   */
  axom_map::Pair<Key, T> insert(Key key, T val)
  {
    axom_map::Bucket<Key, T>* target = &(m_buckets[bucket(get_hash(key))]);
    //Candidate to get cut out if branching becomes too much of an issue.
    axom_map::Pair<Key, T> ret = target->insert_no_update(key, val);
    if(ret.second == true)
    {
      m_size++;
    }
    return ret;
  }
  
  /*!
   * \brief Finds reference to value of key-value pair mapped to supplied pair, and performs insertion if
   *  item does not already exist in Map instance. For reading, modification, and insertion into Map instance.
   *
   * \param [in] key the key used to index into the Map to find/store item
   *
   * \note Major deviation from STL unordered_map in that this function is unsafe, due to inability to internally
   *  call rehash function. User must be careful to ensure integrity of Map with repeated use of this function.
   *
   * \return A reference to the value of key-value pair in the Map instance that is associated with supplied Key.
   *
   */  
  T& operator[](Key key)
  {
    //Since we can't throw an exception, and this needs to return a reference to a variable of type T, 
    //the case where there's no room in the Map is woefully unhandled. 
    axom_map::Pair<Key, T> ins_result = insert(key, T{0});

    return &(ins_result.first->value);
  }
  
  /*!
   * \brief Removes key-value pair with given key from the Map instance.
   *
   * \param [in] key the key of the item to be removed from the Map instance.
   *
   * \note Deviation from STL unordered_map in returning a boolean instead of an iterator to next item in map.
   *  Mostly useful on sparing logic when this moves to paralellism, since the bucket boundaries would become inconvenient.
   *
   * \return A bool value of true if the item was successfully removed, false otherwise. 
   */
  bool erase(Key key)
  {
    axom_map::Bucket<Key, T>* target = &(m_buckets[bucket(get_hash(key))]);
    //Candidate to get cut out if branching becomes too much of an issue.
    bool ret = target->remove(key);
    if(ret == true)
    {
      m_size--;
    }
    return ret;
  }

  /*!
   * \brief Returns key-value pair with given key from the Map instance.
   *
   * \param [in] key the key of the item to be searched for in the Map instance.
   *
   * \return A reference to the requested item if found, sentinel node end otherwise.
   */
  axom_map::Node<Key, T>& find(Key key)
  {
    axom_map::Bucket<Key, T>* target = &(m_buckets[bucket(get_hash(key))]);
    return target->find(key);
  }

  /*!
   * \brief Returns id of bucket associated with a 64-bit integer, intended to be the hash 
   *  of a key.
   *
   * \param [in] hash the id of the bucket to be queried.
   *
   * \note Likely better off inlined, included here in case of later changes to hash collision resolution method.
   *
   * \return A pointer to the queried bucket.
   */
  std::size_t bucket(std::size_t hash) { return hash % m_bucket_count; }

  /*!
   * \brief Returns the maximum number of items per bucket.
   * \return bucket_len the maximum number of items per bucket.
   */
  int bucket_size() { return m_bucket_len; }
  /*!
   * \brief Returns the number of buckets in the Map.
   * \return bucket_count 
   */
  int bucket_count() { return m_bucket_count; }
  /*!
   * \brief Returns the amount of items in the Map instance.
   * \return size the amount of items in the Map instance.
   */
  int size() { return m_size; }
  /*!
   * \brief Returns the overall capacity of the Map instance.
   * \return capacity the overall capacity of the Map instance.
   */
  int max_size() { return m_bucket_len * m_bucket_count; }
  ///@}
private:
  /// \name Private Map Methods
  ///@{
  /*!
   * \brief Memory allocation handler, for Map construction, both at initialization and rehash. Allocates and
   *  initializes linked lists used for data storage, and returns pointer to said array, to be used elsewhere.
   *  Kept separate from constructor and rehasher for sake of modularity in case this needs to be swapped
   *  with a different method in the future.
   *
   * \param [in] bucount the number of buckets to be in side the allocated Map instance. 
   * \param [in] bucklen the amount of items that can fit in a bucket in the allocated Map instance.
   * 
   * \return A pointer to the now-allocated array of linked lists.
   */
  axom_map::Bucket<Key, T>* alloc_map(int bucount, int bucklen)
  {
    axom_map::Bucket<Key, T>* tmp =
      axom::allocate<axom_map::Bucket<Key, T>>(bucount);
    for(int i = 0; i < bucount; i++)
    {
      tmp[i].init(bucklen);
    }

    //update when we're testing our returns
    return tmp;
  }

  /*!
   * \brief Returns hash value for a given input. 
   *
   * \note Currently uses std::hash, will switch to a custom hasher in future.
   *
   * \return A 64-bit integer containing the hashed value of the key.
   */
  std::size_t get_hash(Key input)
  {
    std::size_t hashed = std::hash<Key> {}(input);

    return hashed;
  }

  /// @}

  /// \name Private Data Members
  /// @{

  axom_map::Bucket<Key, T>* m_buckets; /*!< array of pointers to linked lists containing data */
  std::size_t m_bucket_count; /*!< the number of buckets in the Map instance */
  int m_bucket_len; /*!< the number of items that can be contained in a bucket in this Map instance */
  int m_size; /*!< the number of items currenty stored in this Map instance */
  int m_load_factor; /*!< currently unused value, used in STL unordered_map to determine when to resize, which we don't do internally at the moment */

  /// @}
};
} /* namespace experimental */
} /* namespace axom */

#endif /* AXOM_MAP_HPP_ */
