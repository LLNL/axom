/**
 * \file NullSet.h
 *
 * \brief Basic API for a set of entities in a simulation
 *
 */

#ifndef MESHAPI_NULL_SET_H_
#define MESHAPI_NULL_SET_H_

#include "common/Utilities.hpp"
#include "meshapi/Set.hpp"

namespace asctoolkit{
namespace meshapi{


  /**
   * \class NullSet
   *
   * \brief An indexed set (a tuple) of entities in a simulation
   */
    class NullSet : public Set
    {
    public:
        NullSet() {}

        inline SizeType size() const { return SizeType();}

        inline SetIndex at(SetPosition pos) const { verifyPosition(pos); return SetPosition();}
        inline SetIndex operator[](SetPosition pos) const { return at(pos);}

        inline bool isSubset() const { return false; }
        const Set* parentSet() const { return this; }

        bool isValid(bool) const { return true;}

        // TODO: Do we need to add iterator stubs here to satisfy some interface?
        //       The result will be invalid, but it may be useful to get the code to compile, or avoid special logic in the code...
        // iterator begin();
        // iterator end();
        // iterator_pair range();

    private:
        void verifyPosition(SetPosition pos) const
        {
            ATK_ASSERT_MSG(false,"Subscripting on NullSet is never valid."
                           << "\n\tAttempted to access item at index " << pos <<".");
        }
    };


#if 0
    /**
     * \brief NullSets are always equal
     * \note Two sets of different types are (currently) considered to be unequal
     */
    inline bool operator==(NullSet const&, NullSet const&) { ATK_WARNING("operator==(NullSet,NullSet)");return true; }
    /**
     * \brief NullSets are always equal
     * \note Two sets of different types are (currently) considered to be unequal
     */
    inline bool operator!=(NullSet const&, NullSet const&) { ATK_WARNING("operator!=(NullSet,NullSet)");return false; }
#endif


} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_SET_H_
