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
    class NullSet : public Set{
    public:
        NullSet() {}
        inline size_type size() const { return size_type();}
        inline SetIndex operator[](SetIndex idx) const {
            ATK_ASSERT_MSG(false,"Subscripting on NullSet is never valid."
                           << "\n\tAttempted to access item at index " << idx <<".");
            return SetIndex();
        }
    };


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


} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_SET_H_
