/**
 * \file StaticConstantRelation.hpp
 *
 * \brief API for a topological relation between two sets in which entities from the first set
 *        can be related to a constant number of entities from the second set
 *
 */

#ifndef MESHAPI_STATIC_CONSTANT_RELATION_HPP_
#define MESHAPI_STATIC_CONSTANT_RELATION_HPP_

#ifndef MESHAPI_STATIC_CONSTANT_RELATION_ITERATOR_USE_PROXY
//  #define MESHAPI_STATIC_CONSTANT_RELATION_ITERATOR_USE_PROXY
#endif


#include <vector>

//#include <iostream>

#include "slic/slic.hpp"
#include "meshapi/OrderedSet.hpp"
#include "meshapi/Relation.hpp"


namespace asctoolkit {
namespace meshapi    {

  class StaticConstantRelation : public Relation
  {
#ifdef MESHAPI_STATIC_CONSTANT_RELATION_ITERATOR_USE_PROXY
  private:
    /**
     * A small helper class to allow double subscripting on the relation
     */
    class SubscriptProxy {
    public:
      SubscriptProxy(RelationVecConstIterator it, SetPosition stride) : m_iter(it), m_stride(stride) {}
      SetPosition const& operator[](SetPosition index) const
      {
        SLIC_ASSERT_MSG( index < m_stride, "Inner array access out of bounds."
            << "\n\tPresented value: " << index
            << "\n\tMax allowed value: " << static_cast<int>(m_stride - 1));
        return m_iter[index];
      }
    private:
      RelationVecConstIterator m_iter;
      SetPosition m_stride;
    };
#endif

  public:

    typedef Relation::SetPosition                                         SetPosition;

    typedef std::vector<SetPosition>                                      RelationVec;
    typedef RelationVec::iterator                                         RelationVecIterator;
    typedef std::pair<RelationVecIterator,RelationVecIterator>            RelationVecIteratorPair;

    typedef RelationVec::const_iterator                                   RelationVecConstIterator;
    typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>  RelationVecConstIteratorPair;


    typedef OrderedSet< policies::RuntimeSizeHolder<Set::PositionType>      // TODO: change this to a compile time size if/when parent is compile time
                      , policies::RuntimeOffsetHolder<Set::PositionType>
                      , policies::StrideOne<Set::PositionType>
                      , policies::STLVectorIndirection<Set::PositionType, Set::ElementType> > RelationSet;

  public:
    StaticConstantRelation (Set* fromSet = &s_nullSet, Set* toSet = &s_nullSet);
    virtual ~StaticConstantRelation(){}
    /**
     * \note TODO: swap this out for data in the datastore
     */
    void                      bindRelationData(const RelationVec & toOffsets, const SetPosition stride = 0);

    RelationVecConstIterator  begin(SetPosition fromSetIndex)       const
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

#ifdef MESHAPI_STATIC_CONSTANT_RELATION_ITERATOR_USE_PROXY
    const SubscriptProxy operator[](SetPosition fromSetElt) const
    {
      return SubscriptProxy( begin(fromSetElt), size(fromSetElt) );
    }
#else
    /**
     * This function returns the OrderedSet of all elements in the toSet related to 'fromSetElt' in the fromSet.
     */
    const RelationSet operator[](SetPosition fromSetElt) const
    {
        // Note -- we need a better way to initialize an indirection set
        RelationSet rel(size(fromSetElt),toSetBeginIndex(fromSetElt) );
        rel.data() = &m_toSetIndicesVec;

        return rel;
    }


#endif

    SetPosition size(SetPosition fromSetIndex = 0)                  const
    {
      verifyPosition(fromSetIndex);
      return m_stride;
    }

    bool isValid(bool verboseOutput = false) const;

  public:

    /**
     * \name DirectDataAccess
     * \brief Accessor functions to get the underlying relation data
     * \note We will have to figure out a good way to limit this access to situations where it makes sense.
     */

    /// \{
    /**
     * \brief Helper function to access the underlying relation data
     * \note The relation currently 'owns' the underlying vector.
     *       This will be changing soon, and we will only have a reference/pointer to the data.
     */
    RelationVec &       toSetPositionsData()       { return m_toSetIndicesVec; }

    /**
     * \brief Helper function to access the underlying relation data
     * \note The relation currently 'owns' the underlying vector.
     *       This will be changing soon, and we will only have a reference/pointer to the data.
     */
    const RelationVec & toSetPositionsData() const { return m_toSetIndicesVec; }

    /// \}
  private:
    inline void         verifyPosition(SetPosition fromSetIndex)    const { SLIC_ASSERT( fromSetIndex < m_fromSet->size() ); }
    inline SetPosition  toSetBeginIndex(SetPosition fromSetIndex)   const { return m_stride * (fromSetIndex); }
    inline SetPosition  toSetEndIndex(SetPosition fromSetIndex)     const { return m_stride * (fromSetIndex + 1); }



  private:

    SetPosition m_stride;

    Set* m_fromSet;
    Set* m_toSet;

    RelationVec m_toSetIndicesVec;            // vector of toSet entries
  };


} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_STATIC_CONSTANT_RELATION_HPP_
