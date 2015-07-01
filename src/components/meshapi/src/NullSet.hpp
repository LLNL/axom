/**
 * \file NullSet.h
 *
 * \brief Sentinel set type indicating an empty set.
 *
 */

#ifndef MESHAPI_NULL_SET_H_
#define MESHAPI_NULL_SET_H_

#include "slic/slic.hpp"
#include "meshapi/Set.hpp"

namespace asctoolkit {
namespace meshapi {


/**
 * \class NullSet
 *
 * \brief An indexed set (a tuple) of entities in a simulation
 */
  class NullSet : public Set
  {
  public:
    NullSet() {}

    inline PositionType         size() const { return PositionType(); }

    inline ElementType          at(PositionType pos) const { verifyPosition(pos); return PositionType(); }
    inline ElementType operator [](PositionType pos) const { return at(pos); }

    inline bool                 isSubset() const { return false; }
    const Set*                  parentSet() const { return this; }

    bool                        isValid(bool) const { return true; }

    bool                        isEmpty() const { return true; }

    // TODO: Do we need to add iterator stubs here to satisfy some interface?
    //       The result will be invalid, but it may be useful to get the code to compile, or avoid special logic in the code...
    // iterator begin();
    // iterator end();
    // iterator_pair range();

  private:
    void verifyPosition(PositionType pos) const
    {
      SLIC_ASSERT_MSG(false,"Subscripting on NullSet is never valid."
          << "\n\tAttempted to access item at index " << pos << ".");
    }
  };


#if 0
/**
 * \brief NullSets are always equal
 * \note Two sets of different types are (currently) considered to be unequal
 */
  inline bool operator==(NullSet const&, NullSet const&) { ATK_WARNING("operator==(NullSet,NullSet)"); return true; }
/**
 * \brief NullSets are always equal
 * \note Two sets of different types are (currently) considered to be unequal
 */
  inline bool operator!=(NullSet const&, NullSet const&) { ATK_WARNING("operator!=(NullSet,NullSet)"); return false; }
#endif


} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_SET_H_
