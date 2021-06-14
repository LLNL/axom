// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MAP_HPP_
#define AXOM_MAP_HPP_

#include "axom/config.hpp"                 // for compile-time defines
#include "axom/core/Macros.hpp"            // for axom macros
#include "axom/core/memory_management.hpp" // for memory allocation functions
#include "axom/core/Types.hpp"             // for axom types


// C/C++ includes
#include <functional> //for hashing until a custom method is defined
#include <iostream>

namespace axom
{

namespace axom_map
{ 

  /// \name Map Supporting Data Structures
  /// @{
  /*!
   * \brief Node is the foundational data type of the Map, being the lowest-level data storage unit,
   *  with the key and value data for a given item, along with a next value for the linked list.
   *
   */  
  template <typename Key, typename T>
  struct Node{
      IndexType next;
      Key key;
      T value;
      bool operator==(const Node& rhs){
        return (this->next == rhs.next) && (this->key == rhs.key) && (this->value == rhs.value);
      }
            
    }; 

  
  template <typename Key, typename T>
  struct Pair{
    axom_map::Node<Key, T> * first;
    bool second;

    Pair(axom_map::Node<Key, T> * node, bool status) : first(node), second(status) {}
  };
 
  template <typename Key, typename T>   
  class Bucket{
  public:

    Bucket(){}    

    Bucket(int len){
      m_list = axom::allocate <axom_map::Node<Key, T> > (len);
      for(int i = 0; i < len-1; i++){
        m_list[i].next = i+1;
      }
      m_list[len-1].next = -1;
      m_free = 0;
      m_head = -1;
      m_size = 0;
      m_capacity = len;
      m_end.key = Key{0};
      m_end.value = T{0};
      m_end.next = -2;
    }
   
    void init(int len){
      m_list = axom::allocate <axom_map::Node<Key, T> > (len);
      for(int i = 0; i < len-1; i++){
        m_list[i].next = i+1;
      }
      m_list[len-1].next = -1;
      m_free = 0;
      m_head = -1; 
      m_capacity = len;
      m_size = 0;
      m_end.key = Key{0};
      m_end.value = T{0};
      m_end.next = -2;
    }

    axom_map::Pair<Key, T> insert_no_update(Key key, T value){
      if(m_free != -1){
        IndexType ind = m_head;
        if(m_head == -1){
          m_head = m_free;
          ind = m_free;
          m_free = m_list[m_free].next;
          m_list[m_head].next = -1;
          m_list[m_head].key = key;
          m_list[m_head].value = value;
          m_size++;
        }
        else{
          while(m_list[ind].next != -1){
            if(m_list[ind].key == key){
              //change to third value in final version
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
  
   /*bool insert_update(Key key, T value){
     if(free != -1){
        if(head == -1){
          head = free;
          free = list[free].next;
          list[head].next = -1;
          list[head].key = key;
          list[head].value = value;
        }
        else{
          IndexType ind = head;
          while(list[ind].next != -1){
            if(list[ind].key == key){
              //change to third value in final version
              list[ind].value = value;
              return true;
            }
            ind = list[ind].next;
          }
      
          list[ind].next = free;
          ind = free;
          free = list[free].next;
          list[ind].next = -1;
          list[ind].key = key;
          list[ind].value = value;
          size++;
        }
        return true;
      }
      return false;

   }*/
  
    bool remove(Key key){
     if(m_head == -1){
       //add third state for final impl
       return false;
     }   
     IndexType ind = m_head;
     IndexType prev = -1;
     do{
       if(m_list[ind].key == key){
         if(prev != -1){
           m_list[prev].next = m_list[ind].next;
         }
         if(ind == m_head){m_head = m_list[ind].next;}
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

        


   //Pair may be better than a status flag argument. 
   //This is very C.
   axom_map::Node<Key, T>& find(Key key){
     if(m_head == -1){
       return m_end;
     }
     IndexType ind = m_head;
     do{
       if(m_list[ind].key == key){
         return m_list[ind];
       }

       ind = m_list[ind].next;
    }while(ind != -1);
    return m_end;
   }
   
   void print_all(){
     if(m_head == -1){
       return;
     }
     IndexType ind = m_head;
     do{
       std::cout << m_list[ind].value << std::endl;
       ind = m_list[ind].next;
       std::cout << "next step is " << ind << std::endl;
     }while(ind != -1);
   }   
  
  int get_capacity(){ return m_capacity; }
  int get_size(){ return m_size; }   
       
    axom_map::Node<Key, T> * m_list;   
    axom_map::Node<Key, T> m_end;
    IndexType m_head;
    IndexType m_free;
    int m_capacity, m_size;
  
  }; 

  /// @}
}

/*!
 * \class Map
 *
 * \brief Provides a hashmap implementation, relying upon chaining.
 *  Simple hashmap for now, future work will use RAJA and Umpire
 *  to allow use of map within kernels. Resembles unordered_map, 
 *  which we can't use due to a lack of host-device decoration.
 *
 * \tparam Key the type of keys.
 * \tparam T the type of values to hold.
 * \tparam Hash functor that takes an object of Key type and returns 
 *  a hashed value of size_t (for now)
 */ 
 
template <typename Key, typename T, typename Hash = std::hash<Key> >
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
  Map(int num_buckets, int bucket_len=10 ){
    m_bucket_count = buckets;
    m_bucket_len = bucket_len;
    m_buckets = alloc_map(m_bucket_count, m_bucket_len);    
  }

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
  void rehash(int buckets = -1, int factor = -1){
   
    int newlen = 0;
    IndexType ind;
    std::size_t hashed;
    if(buckets != -1 && buckets > m_bucket_count){
      newlen = buckets;
    }
    else if(factor != -1 && factor > 1){
      newlen = factor*m_bucket_count;
    }
    else{
      newlen = 2*m_bucket_count;
    }

    axom_map::Bucket<Key,T> * new_list = alloc_map(newlen, m_bucket_len);

    for(int i = 0; i < m_bucket_count; i++){
      ind = m_buckets[i].m_head;
      while(ind != -1){
        hashed = (get_hash(m_buckets[i].m_list[ind].key) % newlen);
        new_list[hashed].insert_no_update(m_buckets[i].m_list[ind].key, m_buckets[i].m_list[ind].value);     
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
  axom_map::Pair<Key, T> insert(Key key, T val){
    axom_map::Bucket<Key,T> * target = get_bucket(get_hash(key)); 
    //Candidate to get cut out if branching becomes too much of an issue.
    axom_map::Pair<Key, T> ret = target->insert_no_update(key, val);
    if(ret->second == true){ m_size++; }
    return ret;
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
  bool remove(Key key){
    axom_map::Bucket<Key,T> * target = get_bucket(get_hash(key));
    //Candidate to get cut out if branching becomes too much of an issue.
    bool ret = target->remove(key);
    if(ret == true){ m_size--;}
    return ret;
  } 

  /*!
   * \brief Returns key-value pair with given key from the Map instance.
   *
   * \param [in] key the key of the item to be searched for in the Map instance.
   *
   * \return A reference to the requested item if found, sentinel node end otherwise.
   */  
  axom_map::Node<Key, T>& find(Key key){
    axom_map::Bucket<Key,T> * target = get_bucket(get_hash(key));
    return target->find(key);
  }

  /*!
   * \brief Returns the maximum number of items per bucket.
   * \return bucket_len the maximum number of items per bucket.
   */ 
  int get_bucket_size(){ return m_bucket_len; }
  /*!
   * \brief Returns the number of buckets in the Map.
   * \return bucket_count 
   */  
  int get_bucket_count(){ return m_bucket_count; }
  /*!
   * \brief Returns the amount of items in the Map instance.
   * \return size the amount of items in the Map instance.
   */   
  int get_size() { return m_size; }
  /*!
   * \brief Returns the overall capacity of the Map instance.
   * \return capacity the overall capacity of the Map instance.
   */  
  int get_capacity() { return m_bucket_len*m_bucket_count; }

private:
 
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
  axom_map::Bucket<Key,T> * alloc_map(int bucount, int bucklen){
    axom_map::Bucket<Key,T> *tmp = axom::allocate <axom_map::Bucket<Key,T> > (bucount); 
    for(int i = 0; i < bucount; i++){
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
  std::size_t get_hash(Key input){
    std::size_t hashed =  std::hash<Key>{}(input);
     
    return hashed;
  }

  /*!
   * \brief Returns pointer to bucket associated with a 64-bit integer, intended to be the hash 
   *  of a key.
   *
   * \param [in] hash the id of the bucket to be queried.
   *
   * \note Likely better off inlined, included here in case of later changes to hash collision resolution method.
   *
   * \return A pointer to the queried bucket.
   */  
  axom_map::Bucket<Key,T> * get_bucket(std::size_t hash){
    return &(m_buckets[hash%m_bucket_count]);
  }

  /// @}
  
  /// \name Private Data Members
  /// @{
  
  axom_map::Bucket<Key,T> *m_buckets; /*!< array of pointers to linked lists containing data */
  int m_bucket_count; /*!< the number of buckets in the Map instance */
  int m_bucket_len; /*!< the number of items that can be contained in a bucket in this Map instance */
  int m_size; /*!< the number of items currenty stored in this Map instance */
  int m_load_factor; /*!< currently unused value, used in STL unordered_map to determine when to resize, which we don't do internally at the moment */

  /// @}
 

};

} /* namespace axom */

#endif /* AXOM_MAP_HPP_ */
