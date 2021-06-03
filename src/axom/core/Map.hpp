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
  template <typename Key, typename T>
  struct node{
      IndexType next;
      Key key;
      T value;
      bool operator==(const node& rhs){
        return (this->next == rhs.next) && (this->key == rhs.key) && (this->value == rhs.value);
      }
            
    }; 
 
  template <typename Key, typename T>   
  class bucket{
  public:

    bucket(){}    

    bucket(int len){
      list = axom::allocate <axom_map::node<Key, T> > (len);
      for(int i = 0; i < len-1; i++){
        list[i].next = i+1;
      }
      list[len-1].next = -1;
      free = 0;
      head = -1;
      end.key = Key{0};
      end.value = T{0};
      end.next = -2;
    }
   
    void init(int len){
      list = axom::allocate <axom_map::node<Key, T> > (len);
      for(int i = 0; i < len-1; i++){
        list[i].next = i+1;
      }
      list[len-1].next = -1;
      free = 0;
      head = -1;
      end.key = Key{0};
      end.value = T{0};
      end.next = -2;
    }

    bool insert_noupdate(Key key, T value){
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
              return false;
            }
            ind = list[ind].next;
          }
          list[ind].next = free;
          ind = free;
          free = list[free].next;
          list[ind].next = -1;
          list[ind].key = key;
          list[ind].value = value;
        }
        return true;
      }
      return false;
    }
  
   bool insert_update(Key key, T value){
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
        }
        return true;
      }
      return false;

   }
  
   bool remove(Key key){
     if(head == -1){
       //add third state for final impl
       return false;
     }   
     IndexType ind = head;
     IndexType prev = -1;
     do{
       if(list[ind].key == key){
         if(prev != -1){
           list[prev].next = list[ind].next;
         }
         if(ind == head){head = list[ind].next;}
         list[ind].next = free;
         free = ind;
         return true;
       }
       prev = ind;
       ind = list[ind].next;
     } while(ind != -1);
     return false;
   }

        


   //Pair may be better than a status flag argument. 
   //This is very C.
   axom_map::node<Key, T>& find(Key key){
     if(head == -1){
       return end;
     }
     IndexType ind = head;
     do{
       if(list[ind].key == key){
         return list[ind];
       }

       ind = list[ind].next;
    }while(ind != -1);
    return end;
   }
   
   void printall(){
     if(head == -1){
       return;
     }
     IndexType ind = head;
     do{
       std::cout << list[ind].value << std::endl;
       ind = list[ind].next;
       std::cout << "next step is " << ind << std::endl;
     }while(ind != -1);
   }   
  
    
    
    axom_map::node<Key, T> * list;   
    axom_map::node<Key, T> end;
    IndexType head;
    IndexType free;
  
  }; 
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
  Map(int buckets, int bucklen ){
    len = buckets;
    alloc_map(len, bucklen);    
   
   
   insert(27, 236);
   insert(2, 259);
   insert(3, 10);
   insert(35, 2195);
   std::cout << find(3).value;
   remove(3);
   if(find(3) == testlist[0].end){
     std::cout << "removed" << std::endl;
   }
   if(find(4) == testlist[0].end){
     std::cout << "never added" << std::endl;
   }
   std::cout << find(2).value;
  }

  

private:
  
  bool alloc_map(int bucount, int buclen){
    testlist = axom::allocate <axom_map::bucket<Key,T> > (len); 
    for(int i = 0; i < len; i++){
     testlist[i].init(bucklen);
    }

    //update when we're testing our returns
    return true;
  }
   
  axom_map::bucket<Key,T> *testlist;
  int len;
 
  std::size_t get_hash(Key input){
    std::size_t hashed =  std::hash<Key>{}(input);
     
    return hashed;
  }

  axom_map::bucket<Key,T> * get_bucket(std::size_t hash){
    return &(testlist[hash%len]);
  }

public:

  bool insert(Key key, T val){
    axom_map::bucket<Key,T> * target = get_bucket(get_hash(key));  
    return target->insert_noupdate(key, val);
  }

  bool remove(Key key){
    axom_map::bucket<Key,T> * target = get_bucket(get_hash(key));
    return target->remove(key);
  } 

   
  axom_map::node<Key, T>& find(Key key){
    axom_map::bucket<Key,T> * target = get_bucket(get_hash(key));
    return target->find(key);
  }
};

} /* namespace axom */

#endif /* AXOM_MAP_HPP_ */
