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
      size = 0;
      capacity = len;
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
      capacity = len;
      size = 0;
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
          size += 1;
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
          size++;
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
          size++;
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
  
  int get_capacity(){ return capacity; }
  int get_size(){ return size; }   
       
    axom_map::node<Key, T> * list;   
    axom_map::node<Key, T> end;
    IndexType head;
    IndexType free;
    int capacity, size;
  
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
  Map(int buckets, int bucklen=10 ){
    bucket_count = buckets;
    bucket_len = bucklen;
    testlist = alloc_map(bucket_count, bucklen);    
    
   
   insert(27, 236);
   insert(2, 259);
   insert(3, 10);
   insert(35, 2195);
   std::cout << find(3).value << std::endl;
   remove(3);
   if(find(3) == testlist[0].end){
     std::cout << "removed" << std::endl;
   }
   if(find(4) == testlist[0].end){
     std::cout << "never added" << std::endl;
   }
   std::cout << find(2).value << std::endl;
   rehash();
   insert(3,10);
   std::cout << find(3).value << std::endl;
   remove(3);
   if(find(3) == testlist[0].end){
     std::cout << "removed" << std::endl;
   }
   if(find(4) == testlist[0].end){
     std::cout << "never added" << std::endl;
   }
   std::cout << find(2).value << std::endl;


  }

  void rehash(int buckets = -1, int factor = -1){
   
    int newlen = 0;
    IndexType ind;
    std::size_t hashed;
    if(buckets != -1){
      newlen = buckets;
    }
    else if(factor != -1){
      newlen = factor*bucket_count;
    }
    else{
      newlen = 2*bucket_count;
    }

    axom_map::bucket<Key,T> * new_list = alloc_map(newlen, bucket_len);

    for(int i = 0; i < bucket_count; i++){
      ind = testlist[i].head;
      while(ind != -1){
        hashed = (get_hash(testlist[i].list[ind].key) % newlen);
        new_list[hashed].insert_noupdate(testlist[i].list[ind].key, testlist[i].list[ind].value);     
        ind = testlist[i].list[ind].next;
      }
      axom::deallocate(testlist[i].list);
    }
    axom::deallocate(testlist);
    testlist = new_list;
    bucket_count = newlen;
    return; 
  } 

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

  int get_bucket_size(){ return bucket_len; }
  int get_bucket_count(){ return bucket_count; }
  int get_size() { return size; }
  int get_capacity() { return bucket_len*bucket_count; }

private:
  
  axom_map::bucket<Key,T> * alloc_map(int bucount, int bucklen){
    axom_map::bucket<Key,T> *tmp = axom::allocate <axom_map::bucket<Key,T> > (bucount); 
    for(int i = 0; i < bucount; i++){
     tmp[i].init(bucklen);
    }

    //update when we're testing our returns
    return tmp;
  }
   
   
  std::size_t get_hash(Key input){
    std::size_t hashed =  std::hash<Key>{}(input);
     
    return hashed;
  }

  axom_map::bucket<Key,T> * get_bucket(std::size_t hash){
    return &(testlist[hash%bucket_count]);
  }

  axom_map::bucket<Key,T> *testlist;
  int bucket_count, bucket_len, size, load_factor;
 

};

} /* namespace axom */

#endif /* AXOM_MAP_HPP_ */
