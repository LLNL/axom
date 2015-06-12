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
            SubscriptProxy(RelationVecConstIterator it, SetPosition stride): m_iter(it), m_stride(stride) {}
            SetPosition const& operator[](SetPosition index) const
            {
                ATK_ASSERT_MSG( index < m_stride, "Inner array access out of bounds."
                                             <<"\n\tPresented value: "<< index
                                             <<"\n\tMax allowed value: " << static_cast<int>(m_stride -1));
                return m_iter[index];
            }
        private:
            RelationVecConstIterator m_iter;
            SetPosition m_stride;
        };

    public:

        typedef Relation::SetPosition                                           SetPosition;

        typedef std::vector<SetPosition>                                           RelationVec;
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
        void bindRelationData(const RelationVec & toOffsets, const SetPosition stride = 0);

        RelationVecConstIterator begin(SetPosition fromSetIndex)       const
        {
            verifyPosition(fromSetIndex);
            return m_toSetIndicesVec.begin() + toSetBeginIndex(fromSetIndex);
        }

        RelationVecConstIterator end(SetPosition fromSetIndex)         const
        {
            verifyPosition(fromSetIndex);
            return m_toSetIndicesVec.begin() + toSetEndIndex(fromSetIndex);
        }

        RelationVecConstIteratorPair range(SetPosition fromSetIndex)   const
        {
            return std::make_pair(begin(fromSetIndex), end(fromSetIndex));
        }

        SubscriptProxy const operator[](SetPosition fromSetElt) const
        {
            return SubscriptProxy( begin(fromSetElt), size(fromSetElt) );
        }

        SetPosition size(SetPosition fromSetIndex = 0)                  const
        {
            verifyPosition(fromSetIndex);
            return m_stride;
        }

        bool isValid(bool verboseOutput = false) const;

    private:
        inline void  verifyPosition(SetPosition fromSetIndex)       const { ATK_ASSERT( fromSetIndex < m_fromSet->size() ); }
        inline SetPosition toSetBeginIndex(SetPosition fromSetIndex)   const { return m_stride * (fromSetIndex); }
        inline SetPosition toSetEndIndex(SetPosition fromSetIndex)     const { return m_stride * (fromSetIndex+1); }



    private:

        SetPosition       m_stride;

        Set* m_fromSet;
        Set* m_toSet;

        RelationVec m_toSetIndicesVec;        // vector of toSet entries
    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_STATIC_CONSTANT_RELATION_HPP_
