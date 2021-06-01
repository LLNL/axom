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
  Map(int len){
    testlist = bucket(10); 
    testlist.insert_noupdate(37, 1);
    testlist.insert_update(38, 2);
    testlist.remove(38);
    testlist.insert_update(98, 70);
    testlist.insert_update(2, 52);
    testlist.insert_update(24, 12);
    testlist.remove(37);
    testlist.printall();
  }


private:
  
  
  class bucket{
  public:

    bucket(){}    

    bucket(int len){
      list = axom::allocate <node> (len);
      for(int i = 0; i < len-1; i++){
        list[i].next = i+1;
      }
      list[len-1].next = -1;
      free = 0;
      head = -1;
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
   T& find(Key key){
     if(head == -1){
       return list[0].value;
     }
     IndexType ind = head;
     do{
       if(list[ind].key == key){
         return list[ind].value;
       }

       ind = list[ind].next;
    }while(ind != -1);

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
  
    struct node{
      IndexType next;
      Key key;
      T value;
    };     
    
    node * list;   
    node end;
    IndexType head;
    IndexType free;
  
  }; 
  bucket testlist;
};

} /* namespace axom */

#endif /* AXOM_MAP_HPP_ */
