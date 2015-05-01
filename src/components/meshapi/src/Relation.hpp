/**
 * \file Relation.h
 *
 * \brief Basic API for a topological relation between two sets
 *
 */

#ifndef MESHAPI_RELATION_HPP_
#define MESHAPI_RELATION_HPP_

#include <vector>

#include <iostream>

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
/*
            std::cout<<"\nIn size function for fromSetIndex " << fromSetIndex
                    <<"\n\t toSetEndIndex -- "      << toSetEndIndex(fromSetIndex)
                    <<"\n\t toSetBeginIndex -- "    << toSetBeginIndex(fromSetIndex)
                    <<"\n\t size (diff)    -- "     << toSetEndIndex(fromSetIndex) - toSetBeginIndex(fromSetIndex)
                    <<"\n\t fromArray end -- "      << m_fromSetBeginsVec[fromSetIndex+1]
                    <<"\n\t fromArray begin -- "    << m_fromSetBeginsVec[fromSetIndex]
                    <<"\n\t size (diff)    -- "     << m_fromSetBeginsVec[fromSetIndex+1] -m_fromSetBeginsVec[fromSetIndex]
                  <<"\n\t toSetEndIndex cast-- "      << static_cast<unsigned int>(toSetEndIndex(fromSetIndex))
                  <<"\n\t toSetBeginIndex -- "        << static_cast<unsigned int>(toSetBeginIndex(fromSetIndex))
                  <<"\n\t size (diff)    -- "         << static_cast<unsigned int>(toSetEndIndex(fromSetIndex)) - static_cast<unsigned int>(toSetBeginIndex(fromSetIndex))
                    << std::endl;
*/

            return toSetEndIndex(fromSetIndex) - toSetBeginIndex(fromSetIndex);
        }

        bool isValid(bool verboseOutput = false) const;

    private:
        inline void  verifyIndex(Index fromSetIndex)       const { ASSERT( m_fromSet && (fromSetIndex < m_fromSet->size() ) ); }
        inline Index toSetBeginIndex(Index fromSetIndex)   const { return m_fromSetBeginsVec[fromSetIndex]; }
        inline Index toSetEndIndex(Index fromSetIndex)     const { return m_fromSetBeginsVec[fromSetIndex+1]; }



    private:

        OrderedSet* m_fromSet;
        OrderedSet* m_toSet;

        RelationVec m_fromSetBeginsVec;       // vector of size m_fromSet.size() + 1 that points into the to_set vectors
        RelationVec m_toSetIndicesVec;        // vector of toSet entries
    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif /* MESHAPI_RELATION_HPP_ */
