// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file ItemCollection.hpp
 *
 * \brief   Header file for Collection classes.
 *
 *          Each of these classes holds a collection of items of a fixed
 *          type that can be accessed by string name or sidre::IndexType.
 *
 *          The primary intent is to decouple the implementation of the
 *          collections from the Group class which owns collections of
 *          View and child Group objects. They may have other uses,
 *          so they are not dependent on the Group class. Each class is
 *          templated on the item type so that the same class can be used
 *          to hold either View or Group object pointers without
 *          having to code a separate class for each.
 *
 *          By having various collections that obey the same interface,
 *          we can explore alternative collection implementations for
 *          performance (insertion, lookup, etc.) and memory overhead.
 *          The collection used by the Group class can be changed via
 *          the collection type alias in the Group class header file.
 *
 *          To try another collection, encapsulate it in a new class with
 *          the API described below or pass it as a template parameter to
 *          an existing class below if that works.
 *
 *          \attention These classes should be robust against any potential
 *                     user interaction. They don't report errors and leave
 *                     checking of return values to calling code.
 *
 *          \attention Template parameter type must provide a method
 *                     "getName()" that returns a reference to a string object.
 *
 *          \attention The common interface each collection class provides
 *                     is as follows:
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
 *            // (sidre::InvalidName if none).
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
 *          - // Remove item with given name if it exists and return a
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

////////////////////////////////////////////////////////////////////////
//
// ItemCollection keeps an index constant for each item
// as long as it remains in the collection; i.e., don't shift indices
// around.  It has the additional benefit that users can hold on to
// item indices without them being changed without notice.
//
////////////////////////////////////////////////////////////////////////

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
template< typename TYPE >
class ItemCollection
{
public:

  virtual ~ItemCollection() {}

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
