/**
 * \file DynamicVariableRelation.h
 *
 * \brief API for a topological relation between two sets in which entities from the first set
 *        can be related to an arbitrary number of entities from the second set
 *        This relation is dynamic, so it cannot
 *
 */

#ifndef MESHAPI_DYNAMIC_VARIABLE_RELATION_HPP_
#define MESHAPI_DYNAMIC_VARIABLE_RELATION_HPP_

#include <vector>

//#include <iostream>

#include "common/Utilities.hpp"
#include "meshapi/Set.hpp"
#include "meshapi/Relation.hpp"


namespace asctoolkit {
namespace meshapi    {

    class DynamicVariableRelation : public Relation
    {
    public:
        typedef Relation::SetIndex                                              SetIndex;
        typedef Relation::size_type                                             size_type;

        typedef std::vector<SetIndex>                                           RelationVec;
        typedef RelationVec::iterator                                           RelationVecIterator;
        typedef std::pair<RelationVecIterator,RelationVecIterator>              RelationVecIteratorPair;

        typedef RelationVec::const_iterator                                     RelationVecConstIterator;
        typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>    RelationVecConstIteratorPair;

        typedef std::vector< RelationVec>                                       RelationsContainer;
        typedef RelationsContainer::const_iterator                              RelationsContainerCIt;
        typedef RelationsContainer::iterator                                    RelationsContainerIt;

    public:
        DynamicVariableRelation (Set* fromSet = NULL, Set* toSet = NULL);
        virtual ~DynamicVariableRelation(){}

        RelationVecConstIterator begin(SetIndex fromSetIndex)       const
        {
            verifyIndex(fromSetIndex);
            return fromSetRelationsVec(fromSetIndex).begin();
        }

        RelationVecConstIterator end(SetIndex fromSetIndex)         const
        {
            verifyIndex(fromSetIndex);
            return fromSetRelationsVec(fromSetIndex).end();
        }

        RelationVecConstIteratorPair range(SetIndex fromSetIndex)   const
        {
            return std::make_pair(begin(fromSetIndex), end(fromSetIndex));
        }


        RelationVec const& operator[](SetIndex fromSetIndex) const
        {
            verifyIndex(fromSetIndex);
            return m_relationsVec[fromSetIndex];
        }

        size_type size(SetIndex fromSetIndex)                  const
        {
            verifyIndex(fromSetIndex);
            return fromSetRelationsVec(fromSetIndex).size();
        }

        bool isValid(bool verboseOutput = false) const;


    public: // Modifying functions

        void insert(SetIndex fromSetIndex, SetIndex toSetIndex)
        {
            verifyIndex(fromSetIndex);
            m_relationsVec[fromSetIndex].push_back(toSetIndex);

        }

        RelationVec& operator[](SetIndex fromSetIndex)
        {
            verifyIndex(fromSetIndex);
            return m_relationsVec[fromSetIndex];
        }


    private:
        inline void  verifyIndex(SetIndex fromSetIndex)       const { ATK_ASSERT( m_fromSet && (fromSetIndex < m_fromSet->size() ) ); }
        inline RelationVec      & fromSetRelationsVec(SetIndex fromSetIndex)           { return m_relationsVec[fromSetIndex]; }
        inline RelationVec const& fromSetRelationsVec(SetIndex fromSetIndex)   const   { return m_relationsVec[fromSetIndex]; }



    private:

        Set* m_fromSet;
        Set* m_toSet;

        RelationsContainer m_relationsVec;
    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_DYNAMIC_VARIABLE_RELATION_HPP_
