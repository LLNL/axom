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

#include "common/Utilities.hpp"
#include "meshapi/Set.hpp"
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
            SubscriptProxy(RelationVecConstIterator it, SetIndex stride): m_iter(it), m_stride(stride) {}
            SetIndex const& operator[](SetIndex index) const
            {
                ATK_ASSERT_MSG( index < m_stride, "Inner array access out of bounds."
                                             <<"\n\tPresented value: "<< index
                                             <<"\n\tMax allowed value: " << static_cast<int>(m_stride -1));
                return m_iter[index];
            }
        private:
            RelationVecConstIterator m_iter;
            SetIndex m_stride;
        };

    public:

        typedef Relation::SetIndex                                              SetIndex;
        typedef Relation::size_type                                            size_type;

        typedef std::vector<SetIndex>                                           RelationVec;
        typedef RelationVec::iterator                                           RelationVecIterator;
        typedef std::pair<RelationVecIterator,RelationVecIterator>              RelationVecIteratorPair;

        typedef RelationVec::const_iterator                                     RelationVecConstIterator;
        typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>    RelationVecConstIteratorPair;

    public:
        StaticConstantRelation (Set* fromSet = &s_nullSet, Set* toSet = &s_nullSet);
        virtual ~StaticConstantRelation(){}
        /**
         * \note TODO: swap this out for data in the datastore
         */
        void setRelation(RelationVec const& toOffsets, SetIndex stride = 0);

        RelationVecConstIterator begin(SetIndex fromSetIndex)       const
        {
            verifyIndex(fromSetIndex);
            return m_toSetIndicesVec.begin() + toSetBeginIndex(fromSetIndex);
        }

        RelationVecConstIterator end(SetIndex fromSetIndex)         const
        {
            verifyIndex(fromSetIndex);
            return m_toSetIndicesVec.begin() + toSetEndIndex(fromSetIndex);
        }

        RelationVecConstIteratorPair range(SetIndex fromSetIndex)   const
        {
            return std::make_pair(begin(fromSetIndex), end(fromSetIndex));
        }

        SubscriptProxy const operator[](SetIndex fromSetElt) const
        {
            return SubscriptProxy( begin(fromSetElt), size(fromSetElt) );
        }

        size_type size(SetIndex fromSetIndex = 0)                  const
        {
            verifyIndex(fromSetIndex);        
            return m_stride;
        }

        bool isValid(bool verboseOutput = false) const;

    private:
        inline void  verifyIndex(SetIndex fromSetIndex)       const { ATK_ASSERT( fromSetIndex < m_fromSet->size() ); }
        inline SetIndex toSetBeginIndex(SetIndex fromSetIndex)   const { return m_stride * (fromSetIndex); }
        inline SetIndex toSetEndIndex(SetIndex fromSetIndex)     const { return m_stride * (fromSetIndex+1); }



    private:

        SetIndex       m_stride;

        Set* m_fromSet;
        Set* m_toSet;

        RelationVec m_toSetIndicesVec;        // vector of toSet entries
    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_STATIC_CONSTANT_RELATION_HPP_
