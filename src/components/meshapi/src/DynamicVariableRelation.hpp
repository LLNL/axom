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
#include "meshapi/OrderedSet.hpp"
#include "meshapi/Relation.hpp"


namespace asctoolkit {
namespace meshapi    {

    class DynamicVariableRelation : public Relation
    {
    public:
        typedef MeshIndexType                                                   Index;
        typedef MeshSizeType                                                    size_type;

        typedef std::vector<Index>                                              RelationVec;
        typedef RelationVec::iterator                                           RelationVecIterator;
        typedef std::pair<RelationVecIterator,RelationVecIterator>              RelationVecIteratorPair;

        typedef RelationVec::const_iterator                                     RelationVecConstIterator;
        typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>    RelationVecConstIteratorPair;

        typedef std::vector< RelationVec>                                       RelationsContainer;
        typedef RelationsContainer::const_iterator                              RelationsContainerCIt;
        typedef RelationsContainer::iterator                                    RelationsContainerIt;

    public:
        DynamicVariableRelation (OrderedSet* fromSet = NULL, OrderedSet* toSet = NULL);
        virtual ~DynamicVariableRelation(){}

        RelationVecConstIterator begin(Index fromSetIndex)       const
        {
            verifyIndex(fromSetIndex);
            return fromSetRelationsVec(fromSetIndex).begin();
        }

        RelationVecConstIterator end(Index fromSetIndex)         const
        {
            verifyIndex(fromSetIndex);
            return fromSetRelationsVec(fromSetIndex).end();
        }

        RelationVecConstIteratorPair range(Index fromSetIndex)   const
        {
            return std::make_pair(begin(fromSetIndex), end(fromSetIndex));
        }


        RelationVec const& operator[](Index fromSetIndex) const
        {
            verifyIndex(fromSetIndex);
            return m_relationsVec[fromSetIndex];
        }

        size_type size(Index fromSetIndex)                  const
        {
            verifyIndex(fromSetIndex);
            return fromSetRelationsVec(fromSetIndex).size();
        }

        bool isValid(bool verboseOutput = false) const;


    public: // Modifying functions

        void insert(Index fromSetIndex, Index toSetIndex)
        {
            verifyIndex(fromSetIndex);
            m_relationsVec[fromSetIndex].push_back(toSetIndex);

        }

        RelationVec& operator[](Index fromSetIndex)
        {
            verifyIndex(fromSetIndex);
            return m_relationsVec[fromSetIndex];
        }


    private:
        inline void  verifyIndex(Index fromSetIndex)       const { ATK_ASSERT( m_fromSet && (fromSetIndex < m_fromSet->size() ) ); }
        inline RelationVec      & fromSetRelationsVec(Index fromSetIndex)           { return m_relationsVec[fromSetIndex]; }
        inline RelationVec const& fromSetRelationsVec(Index fromSetIndex)   const   { return m_relationsVec[fromSetIndex]; }



    private:

        OrderedSet* m_fromSet;
        OrderedSet* m_toSet;

        RelationsContainer m_relationsVec;
    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_DYNAMIC_VARIABLE_RELATION_HPP_
