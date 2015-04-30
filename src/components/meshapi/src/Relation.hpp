/**
 * \file OrderedSet.h
 *
 * \brief Basic API for a topological relation between two sets
 *
 */

#ifndef MESHAPI_RELATION_HPP_
#define MESHAPI_RELATION_HPP_

#include <vector>
#include "meshapi/Utilities.hpp"
#include "meshapi/OrderedSet.hpp"


namespace asctoolkit {
namespace meshapi    {

    class Relation
    {
    public:
        typedef MeshIndexType                                          Index;
        typedef MeshSizeType                                           size_type;

        typedef std::vector<Index>                                     RelationVec;
        typedef typename RelationVec::iterator                         RelationVecIterator;
        typedef std::pair<RelationVecIterator,RelationVecIterator>     RelationVecIteratorPair;

        typedef typename RelationVec::const_iterator                   RelationVecConstIterator;
        typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>     RelationVecConstIteratorPair;

    public:
        Relation (OrderedSet* fromSet, OrderedSet* toSet);

        /**
         * \note TODO: swap this out for data in the datastore
         */
        void setRelation(RelationVec const& beginsVec, RelationVec const& toOffsets);

        RelationVecConstIterator begin(Index fromSetIndex)       const
        {
            verifyIndex(fromSetIndex);
            return m_toSetIndicesVec.begin() + toSetBeginIndex(fromSetIndex);
        }

        RelationVecConstIterator end(Index fromSetIndex)         const
        {
            verifyIndex(fromSetIndex);
            return m_toSetIndicesVec.begin() + toSetEndIndex(fromSetIndex);
        }

        RelationVecConstIteratorPair range(Index fromSetIndex)   const
        {
            return std::make_pair(begin(fromSetIndex), end(fromSetIndex));
        }


        size_type size(Index fromSetIndex)                  const
        {
            verifyIndex(fromSetIndex);
            return toSetEndIndex(fromSetIndex) - toSetBeginIndex(fromSetIndex);
        }

        bool isValid(bool verboseOutput = false) const;

    private:
        inline void  verifyIndex(Index fromSetIndex)       const { ASSERT( m_fromSet && (fromSetIndex < m_fromSet->size() ) ); }
        inline Index toSetBeginIndex(Index fromSetIndex)   const { m_fromSetBeginsVec[fromSetIndex]; }
        inline Index toSetEndIndex(Index fromSetIndex)     const { m_fromSetBeginsVec[fromSetIndex+1]; }



    private:

        OrderedSet* m_fromSet;
        OrderedSet* m_toSet;

        RelationVec m_fromSetBeginsVec;       // vector of size m_fromSet.size() + 1 that points into the to_set vectors
        RelationVec m_toSetIndicesVec;        // vector of toSet entries
    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif /* MESHAPI_RELATION_HPP_ */
