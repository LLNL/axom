/**
 * \file Relation.h
 *
 * \brief Basic API for a topological relation between two sets
 *
 */

#ifndef MESHAPI_RELATION_HPP_
#define MESHAPI_RELATION_HPP_

#include <vector>


namespace asctoolkit {
namespace meshapi    {

    class Relation
    {
    public:
        typedef MeshIndexType                                          Index;
        typedef MeshSizeType                                           size_type;

        typedef std::vector<Index>                                     RelationVec;
        typedef RelationVec::iterator                         RelationVecIterator;
        typedef std::pair<RelationVecIterator,RelationVecIterator>     RelationVecIteratorPair;

        typedef RelationVec::const_iterator                   RelationVecConstIterator;
        typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>     RelationVecConstIteratorPair;

    public:
        virtual ~Relation(){}

        //void setRelation(RelationVec const& beginsVec, RelationVec const& toOffsets) = 0;

        virtual RelationVecConstIterator begin(Index fromSetIndex)       const  = 0;

        virtual RelationVecConstIterator end(Index fromSetIndex)         const  = 0;

        virtual RelationVecConstIteratorPair range(Index fromSetIndex)   const  = 0;

        virtual size_type size(Index fromSetIndex)                       const  = 0;

        virtual bool isValid(bool verboseOutput = false)                const = 0;

    };


} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_RELATION_HPP_
