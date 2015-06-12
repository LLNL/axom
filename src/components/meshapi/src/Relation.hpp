/**
 * \file Relation.h
 *
 * \brief Basic API for a topological relation between two sets
 *
 */

#ifndef MESHAPI_RELATION_HPP_
#define MESHAPI_RELATION_HPP_

#include <vector>

#include "meshapi/Set.hpp"
#include "meshapi/NullSet.hpp"

namespace asctoolkit {
namespace meshapi    {

    class NullSet;

    class Relation
    {
    public:
        typedef Set::PositionType                                                SetPosition;

        typedef std::vector<SetPosition>                                        RelationVec;
        typedef RelationVec::iterator                                           RelationVecIterator;
        typedef std::pair<RelationVecIterator,RelationVecIterator>              RelationVecIteratorPair;

        typedef RelationVec::const_iterator                                     RelationVecConstIterator;
        typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>    RelationVecConstIteratorPair;

        static NullSet s_nullSet;

    public:
        virtual ~Relation(){}

        //void setRelation(RelationVec const& beginsVec, RelationVec const& toOffsets) = 0;

        virtual RelationVecConstIterator begin(SetPosition fromSetIndex)       const  = 0;

        virtual RelationVecConstIterator end(SetPosition fromSetIndex)         const  = 0;

        virtual RelationVecConstIteratorPair range(SetPosition fromSetIndex)   const  = 0;

        virtual SetPosition size(SetPosition fromSetIndex)                       const  = 0;

        virtual bool isValid(bool verboseOutput = false)                const = 0;

    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_RELATION_HPP_
