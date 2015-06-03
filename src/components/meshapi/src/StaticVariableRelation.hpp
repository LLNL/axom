/**
 * \file StaticVariableRelation.h
 *
 * \brief API for a topological relation between two sets in which entities from the first set
 *        can be related to an arbitrary number of entities from the second set
 *
 */

#ifndef MESHAPI_STATIC_VARIABLE_RELATION_HPP_
#define MESHAPI_STATIC_VARIABLE_RELATION_HPP_

#include <vector>

//#include <iostream>

#include "common/Utilities.hpp"
#include "meshapi/Set.hpp"
#include "meshapi/NullSet.hpp"
#include "meshapi/Relation.hpp"


namespace asctoolkit {
namespace meshapi    {

    class StaticVariableRelation : public Relation
    {
    private:
        /**
         * A small helper class to allow double subscripting on the relation
         */
        class SubscriptProxy{
        public:
            SubscriptProxy(RelationVecConstIterator it, SizeType size): m_iter(it), m_size(size) {}
            SetIndex const& operator[](SetIndex index) const
            {
                ATK_ASSERT_MSG( index < static_cast<SetIndex>(m_size), "Inner array access out of bounds."
                                             <<"\n\tPresented value: "<< index
                                             <<"\n\tMax allowed value: " << static_cast<int>(m_size -1));
                return m_iter[index];
            }
        private:
            RelationVecConstIterator m_iter;
            SizeType m_size;
        };
    public:
        typedef Relation::SetIndex                                      SetIndex;
        typedef Relation::SizeType                                      SizeType;
        typedef Relation::SetPosition                                   SetPosition;

        typedef std::vector<SetIndex>                                   RelationVec;
        typedef RelationVec::iterator                                   RelationVecIterator;
        typedef std::pair<RelationVecIterator,RelationVecIterator>      RelationVecIteratorPair;

        typedef RelationVec::const_iterator                   RelationVecConstIterator;
        typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>     RelationVecConstIteratorPair;



    public:
        StaticVariableRelation (Set* fromSet = &s_nullSet, Set* toSet = &s_nullSet);
        virtual ~StaticVariableRelation(){}
        /**
         * \note TODO: swap this out for data in the datastore
         */
        void setRelation(RelationVec const& beginsVec, RelationVec const& toOffsets);

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

        SizeType size(SetIndex fromSetIndex)                  const
        {
            verifyIndex(fromSetIndex);
            return toSetEndIndex(fromSetIndex) - toSetBeginIndex(fromSetIndex);
        }

        bool isValid(bool verboseOutput = false) const;

    private:
        inline void  verifyIndex(SetIndex fromSetIndex)       const { ATK_ASSERT( fromSetIndex <  static_cast<SetIndex>(m_fromSet->size() ) ); }
        inline SetIndex toSetBeginIndex(SetIndex fromSetIndex)   const { return m_fromSetBeginsVec[fromSetIndex]; }
        inline SetIndex toSetEndIndex(SetIndex fromSetIndex)     const { return m_fromSetBeginsVec[fromSetIndex+1]; }



    private:

        Set* m_fromSet;
        Set* m_toSet;

        RelationVec m_fromSetBeginsVec;       // vector of size m_fromSet.size() + 1 that points into the to_set vectors
        RelationVec m_toSetIndicesVec;        // vector of toSet entries
    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_STATIC_VARIABLE_RELATION_HPP_
