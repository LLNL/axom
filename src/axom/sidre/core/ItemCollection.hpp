// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file ItemCollection.hpp
 *
 * \brief   Header file for ItemCollection.
 *
 *          This is a templated abstract base class defining an interface for
 *          classes holding a collection of items of a fixed
 *          type that can be accessed by string name or sidre::IndexType.
 *
 *          The primary intent is to decouple the implementation of the
 *          collections from the Group class which owns collections of
 *          View and child Group objects. They may have other uses,
 *          so they are not dependent on the Group class. This class is
 *          templated on the item type so that derived classes can be used
 *          to hold either View or Group object pointers without
 *          having to code a separate class for each.
 *
 *          Derived implemenations of this class can be used to explore
 *          alternative collection implementations for performance
 *          (insertion, lookup, etc.) and memory overhead.
 *
 *          \attention These classes should be robust against any potential
 *                     user interaction. They don't report errors and leave
 *                     checking of return values to calling code.
 *
 *          \attention The interface defined by this class is as follows:
 *
 *          \verbatim
 *
 *          - // Return number of items in collection.
 *
 *               size_t getNumItems() const;
 *
 *          - // Return first valid item index for iteration.
 *            // sidre::InvalidIndex returned if no items in collection
 *
 *               IndexType getFirstValidIndex() const;
 *
 *          - // Return next valid item index for iteration.
 *            // sidre::InvalidIndex returned if there are no more items
 *            // to be iterated over.
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

#ifndef SIDRE_ITEMCOLLECTIONS_HPP_
#define SIDRE_ITEMCOLLECTIONS_HPP_

// Other axom headers
#include "axom/config.hpp"
#include "axom/core/Types.hpp"

// Sidre project headers
#include "SidreTypes.hpp"

namespace axom
{
namespace sidre
{
/*!
 *************************************************************************
 *
 * \class ItemCollection
 *
 * \brief ItemCollection is an abstract base class template for holding
 *        a collection of items of template parameter type TYPE.  Derived
 *        child classes can determine how to specifically store the items.
 *
 *************************************************************************
 */
template <typename TYPE>
class ItemCollection
{
public:
  virtual ~ItemCollection() { }

  //
  // Default compiler-generated ctor, dtor, copy ctor, and copy assignment
  // operator suffice for this class.
  //

  ///
  virtual size_t getNumItems() const = 0;

  ///
  virtual IndexType getFirstValidIndex() const = 0;

  ///
  virtual IndexType getNextValidIndex(IndexType idx) const = 0;

  ///
  virtual bool hasItem(const std::string& name) const = 0;

  ///
  virtual bool hasItem(IndexType idx) const = 0;

  ///
  virtual TYPE* getItem(const std::string& name) = 0;

  ///
  virtual TYPE const* getItem(const std::string& name) const = 0;

  ///
  virtual TYPE* getItem(IndexType idx) = 0;

  ///
  virtual TYPE const* getItem(IndexType idx) const = 0;

  ///
  virtual const std::string& getItemName(IndexType idx) const = 0;

  ///
  virtual IndexType getItemIndex(const std::string& name) const = 0;

  ///
  virtual IndexType insertItem(TYPE* item, const std::string& name) = 0;

  ///
  virtual TYPE* removeItem(const std::string& name) = 0;

  ///
  virtual TYPE* removeItem(IndexType idx) = 0;

  ///
  virtual void removeAllItems() = 0;

private:
};

} /* end namespace sidre */
} /* end namespace axom */

#endif /* SIDRE_ITEMCOLLECTIONS_HPP_ */
