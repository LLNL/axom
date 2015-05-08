/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Header file for Collection classes.
 *
 *          They are used to hold a collection of items of a fixed type
 *          that can be accessed by string name or common::IDType index.
 *
 *          These are mainly used by the DataGroup class for holding
 *          collections of DataView and child DataGroup objects. But,
 *          they may have other uses. So they are not dependent on the
 *          DataGroup class.
 *
 *          The primary goal is to have fixed collection interface in the 
 *          DataGroup class implementation and be able to try out alternative
 *          collection implementations for performance (insertion, lookup, 
 *          etc.) and memory overhead.
 *
 *          Each of these classes is a template on the item type so that they 
 *          can be used to hold either DataView or DataGroup object pointers
 *          without having to code a separate class for each.
 *
 *          IMPORTANT: These classes should be robust against any potential
 *                     user interaction. They don't report errors and leave
 *                     checking of return values to calling code.
 *
 *          IMPORTANT: Template parameter type must provide a method 
 *                     "name()" that returns a string object.
 *
 *          IMPORTANT: For API consistency, each collection class must provide
 *                     the following methods:
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
 *               bool hasItem(common::IDType idx) const;
 *
 *          - // Return pointer to item with given name (nullptr if none).
 *
 *               TYPE* getItem(const std::string& name);
 *               TYPE const* getItem(const std::string& name) const ;
 *
 *          - // Return pointer to item with given index (nullptr if none).
 *
 *               TYPE* getItem(common::IDType idx);
 *               TYPE const* getItem(common::IDType idx) cosnt;
 *
 *          - // Insert item with given name; return true if insertion 
 *               succeeded, and false otherwise.
 *
 *               bool insertItem(TYPE* item, const std::string& name);
 *
 *          - // Remove item with given name if it exists and return a 
 *               pointer to it. If it doesn't exist, return nullptr.
 *
 *               TYPE* removeItem(const std::string& name);
 *
 *          - // Remove item with given name if it exists and return a 
 *               pointer to it. If it doesn't exist, return nullptr.
 *
 *               TYPE* removeItem(common::IDType idx);
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
#include "common/Types.hpp"

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
* \warning Only std::map and std::unordered_map have been tested so far.
*        
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

   size_t getNumItems() const 
   { 
      return m_items.size();
   }    

   bool hasItem(const std::string& name) const
   {
      typename MAP_TYPE::const_iterator mit = m_name2idx_map.find(name);
      return ( mit != m_name2idx_map.end() ? true : false );
   }

   bool hasItem(common::IDType idx) const 
   { 
       return (idx < getNumItems() && m_items[idx]); 
   }

   TYPE* getItem(const std::string& name)
   {
      typename MAP_TYPE::iterator mit = m_name2idx_map.find(name);
      return ( mit != m_name2idx_map.end() ? 
               m_items[ mit->second ] : nullptr );
   }

   TYPE const* getItem(const std::string& name) const
   {
      typename MAP_TYPE::const_iterator mit = m_name2idx_map.find(name);
      return ( mit != m_name2idx_map.end() ?
               m_items[ mit->second ] : nullptr );
   }

   TYPE* getItem(common::IDType idx)
   {
      return ( hasItem(idx) ? m_items[idx] : nullptr );
   }

   TYPE const* getItem(common::IDType idx) const
   {
      return ( hasItem(idx) ? m_items[idx] : nullptr );
   }

   bool insertItem(TYPE* item, const std::string& name);

   TYPE* removeItem(const std::string& name);

   TYPE* removeItem(common::IDType idx);

   void removeAllItems()
   {
      m_items.clear();
      m_name2idx_map.clear();
   } 

private:
   std::vector<TYPE*>  m_items;
   MAP_TYPE            m_name2idx_map;
};

template <typename TYPE, typename MAP_TYPE> 
bool MapCollection<TYPE, MAP_TYPE>::insertItem(TYPE* item, 
                                               const std::string& name)
{
   if ( m_name2idx_map.insert( std::make_pair(name, m_items.size()) ).second ) {
      // item was inserted into map
      m_items.push_back( item );
      return true;
   } else {
      // item was NOT inserted into map
      return false;
   }
} 

template <typename TYPE, typename MAP_TYPE>
TYPE* MapCollection<TYPE, MAP_TYPE>::removeItem(const std::string& name)
{
   TYPE* ret_val = nullptr;

   typename MAP_TYPE::iterator mit = m_name2idx_map.find(name);
   if ( mit != m_name2idx_map.end() ) {
      common::IDType idx = mit->second;

      ret_val = m_items[idx];

      m_name2idx_map.erase(mit);
      m_items.erase(m_items.begin() + idx);

      // Decrement approriate item indices
      for (mit = m_name2idx_map.begin(); mit != m_name2idx_map.end(); ++mit)
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
TYPE* MapCollection<TYPE, MAP_TYPE>::removeItem(common::IDType idx)
{
   if ( hasItem(idx) ) {
      TYPE* item = removeItem( m_items[idx]->getName() );
      return item;
   } else {
      return nullptr;
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
