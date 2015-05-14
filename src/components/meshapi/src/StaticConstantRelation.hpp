/**
 * \file StaticConstantRelation.h
 *
 * \brief API for a topological relation between two sets in which entities from the first set
 *        can be related to a constant number of entities from the second set
 *
 */

#ifndef MESHAPI_STATIC_CONSTANT_RELATION_HPP_
#define MESHAPI_STATIC_CONSTANT_RELATION_HPP_

#include <vector>

//#include <iostream>

#include "meshapi/Utilities.hpp"
#include "meshapi/OrderedSet.hpp"
#include "meshapi/Relation.hpp"


namespace asctoolkit {
namespace meshapi    {

    class StaticConstantRelation : public Relation
    {
    private:
        /**
         * A small helper class to allow double subscripting on the relation
         */
        class SubscriptProxy{
        public:
            SubscriptProxy(RelationVecConstIterator it, Index stride): m_iter(it), m_stride(stride) {}
            Index const& operator[](Index index) const
            {
                ASSERT2( index < m_stride, "Inner array access out of bounds."
                                             <<"\n\tPresented value: "<< index
                                             <<"\n\tMax allowed value: " << static_cast<int>(m_stride -1))
                return m_iter[index];
            }
        private:
            RelationVecConstIterator m_iter;
            Index m_stride;
        };

    public:

        typedef MeshIndexType                                                   Index;
        typedef MeshSizeType                                                    size_type;

        typedef std::vector<Index>                                              RelationVec;
        typedef RelationVec::iterator                                           RelationVecIterator;
        typedef std::pair<RelationVecIterator,RelationVecIterator>              RelationVecIteratorPair;

        typedef RelationVec::const_iterator                                     RelationVecConstIterator;
        typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>    RelationVecConstIteratorPair;

    public:
        StaticConstantRelation (OrderedSet* fromSet = NULL, OrderedSet* toSet = NULL);
        virtual ~StaticConstantRelation(){}
        /**
         * \note TODO: swap this out for data in the datastore
         */
        void setRelation(RelationVec const& toOffsets, Index stride = 0);

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

        SubscriptProxy const operator[](Index fromSetElt) const
        {
            return SubscriptProxy( begin(fromSetElt), size(fromSetElt) );
        }

        size_type size(Index fromSetIndex = 0)                  const
        {
            verifyIndex(fromSetIndex);        
            return m_stride;
        }

        bool isValid(bool verboseOutput = false) const;

    private:
        inline void  verifyIndex(Index fromSetIndex)       const { ASSERT( m_fromSet && (fromSetIndex < m_fromSet->size() ) ); }
        inline Index toSetBeginIndex(Index fromSetIndex)   const { return m_stride * (fromSetIndex); }
        inline Index toSetEndIndex(Index fromSetIndex)     const { return m_stride * (fromSetIndex+1); }



    private:

        Index       m_stride;

        OrderedSet* m_fromSet;
        OrderedSet* m_toSet;

        RelationVec m_toSetIndicesVec;        // vector of toSet entries
    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_STATIC_CONSTANT_RELATION_HPP_
