/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/**
 * \file
 *
 * \brief API for a topological relation between two sets in which entities from the first set
 *        can be related to a constant number of entities from the second set
 *
 */

#ifndef SLAM_STATIC_CONSTANT_RELATION_HPP_
#define SLAM_STATIC_CONSTANT_RELATION_HPP_

#ifndef SLAM_STATIC_CONSTANT_RELATION_ITERATOR_USE_PROXY
//  #define SLAM_STATIC_CONSTANT_RELATION_ITERATOR_USE_PROXY
#endif


#include <vector>

//#include <iostream>

#include "common/ATKMacros.hpp"
#include "slic/slic.hpp"

#include "slam/OrderedSet.hpp"
#include "slam/Relation.hpp"
#include "slam/PolicyTraits.hpp"


// TODO: We need to add policies to this relation class
//  We already have:
//  * StridePolicy -- this dictates whether the constant for the striding of each relation will be provided at runtime or compile time
//  We are missing:
//  * The set type of the begins set -- this is determined by the stride policy
//                                   -- more generally, we can combine this set type with the variable static relation by allowing this to be an indirection set (of size fromSet.size()+1
//                                      or a strided set (of size fromSet.size()) without indirection
//
//  * The set type of the offsets set -- this can be an indirection set into ToSet
//                                    -- or we can allow it to be an offset-based position set
//                                      (e.g. for relations of type: Zone-to-Side / Zone-to-Corner -- where rel[i][j] = SZ*i+j -- offset is SZ*i, size of set is SZ, stride is 1, no indirection)
//                                      The storage cost here is practically 0 for the entire relation



namespace asctoolkit {
namespace slam    {

  template< typename StridePolicy = policies::RuntimeStrideHolder<Set::PositionType>
  , typename FromSetType = Set
  , typename ToSetType = Set
  >
  class StaticConstantRelation : public Relation
                                 , StridePolicy
  {
#ifdef SLAM_STATIC_CONSTANT_RELATION_ITERATOR_USE_PROXY
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

    typedef StridePolicy                                                                                              StridePolicyType;

    //-----

    typedef Relation::SetPosition                                                                                     SetPosition;

    typedef std::vector<SetPosition>                                                                                  RelationVec;
    typedef RelationVec::iterator                                                                                     RelationVecIterator;
    typedef std::pair<RelationVecIterator,RelationVecIterator>                                                        RelationVecIteratorPair;

    typedef RelationVec::const_iterator                                                                               RelationVecConstIterator;
    typedef std::pair<RelationVecConstIterator,RelationVecConstIterator>                                              RelationVecConstIteratorPair;

    typedef typename policies::StrideToSize<StridePolicyType, SetPosition, StridePolicyType::DEFAULT_VALUE>::SizeType CorrespondingSizeType;
    typedef OrderedSet< CorrespondingSizeType      // The cardinality of each relational operator is determined by the StridePolicy of the relation
        , policies::RuntimeOffsetHolder<Set::PositionType>
        , policies::StrideOne<Set::PositionType>
        , policies::STLVectorIndirection<Set::PositionType, Set::ElementType> > RelationSet;

    struct RelationBuilder;

  public:
    StaticConstantRelation ( FromSetType* fromSet = EmptySetTraits<FromSetType>::emptySet()
        , ToSetType* toSet = EmptySetTraits<ToSetType>::emptySet() )
        : StridePolicy(CorrespondingSizeType::DEFAULT_VALUE), m_fromSet(fromSet), m_toSet(toSet) {}

    StaticConstantRelation( const RelationBuilder & builder)
        : StridePolicy(builder.m_stride)
          , m_fromSet(builder.m_fromSet)
          , m_toSet(builder.m_toSet)
    {}

    ~StaticConstantRelation(){}


  public:
    struct RelationBuilder
    {
      friend class StaticConstantRelation;

      RelationBuilder()
          : m_fromSet( EmptySetTraits<FromSetType>::emptySet() )
            , m_toSet( EmptySetTraits<ToSetType>::emptySet() )
      {}

      RelationBuilder&  fromSet(FromSetType* pFromSet)  { m_fromSet = pFromSet; return *this; }
      RelationBuilder&  toSet(ToSetType* pToSet)        { m_toSet  = pToSet; return *this; }
      RelationBuilder&  stride(SetPosition str)         { m_stride = StridePolicyType(str); return *this; }

      // This needs to wait until relation data is set by a policy....
      //RelationBuilder& offsets()   {  return *this;}

    private:
      FromSetType* m_fromSet;
      ToSetType* m_toSet;
      StridePolicyType m_stride;
    };


  public:


    /**
     * \note TODO: swap this out for data in the datastore
     */
    void bindRelationData(const RelationVec & toOffsets, const SetPosition stride = StridePolicyType::DEFAULT_VALUE)
    {
      StridePolicyType::setStride(stride);

      m_toSetIndicesVec.clear();
      m_toSetIndicesVec.reserve(toOffsets.size());
      std::copy(toOffsets.begin(), toOffsets.end(), std::back_inserter(m_toSetIndicesVec));
    }

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

#ifdef SLAM_STATIC_CONSTANT_RELATION_ITERATOR_USE_PROXY
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
      typedef typename RelationSet::SetBuilder SetBuilder;
      return SetBuilder()
             .size( elemSize(fromSetElt))
             .offset( toSetBeginIndex(fromSetElt) )
             .data( &m_toSetIndicesVec)
      ;
    }
#endif

    /**
     * Returns the number of elements in toSet related to item at position fromSetIndex of fromSet
     */
    SetPosition size(SetPosition fromSetIndex = 0)                  const
    {
      return elemSize(fromSetIndex);
    }

    /**
     * Checks the validity of the relation
     */
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
    inline SetPosition  elemSize(SetPosition fromSetIndex) const
    {
      verifyPosition(fromSetIndex);
      return stride();
    }
    inline void verifyPosition(SetPosition ATK_DEBUG_PARAM(fromSetIndex))    const
    {
      SLIC_ASSERT_MSG( fromSetIndex >= 0 && fromSetIndex < m_fromSet->size(),
          "Failed verify position with position " << fromSetIndex
                                                  << ". Valid positions are integers between 0 and "
                                                  << (m_fromSet->size() - 1)
      );
    }
    inline SetPosition  toSetBeginIndex(SetPosition fromSetIndex)   const { return stride() * (fromSetIndex); }
    inline SetPosition  toSetEndIndex(SetPosition fromSetIndex)     const { return stride() * (fromSetIndex + 1); }

    inline SetPosition  stride() const { return StridePolicyType::stride(); }

  private:

    Set* m_fromSet;
    Set* m_toSet;

    RelationVec m_toSetIndicesVec;            // vector of toSet entries
  };


  template<typename FromSetType, typename ToSetType, typename StridePolicy>
  bool StaticConstantRelation<FromSetType, ToSetType, StridePolicy>::isValid(bool verboseOutput) const
  {
    bool bValid = true;

    std::stringstream errSstr;

    if( EmptySetTraits<FromSetType>::isEmpty(m_fromSet) || EmptySetTraits<ToSetType>::isEmpty(m_toSet) )
    {
      if(!m_toSetIndicesVec.empty())
      {
        if(verboseOutput)
        {
          errSstr << "\n\t* toSetIndicesVec was not empty "
                  << " -- fromSet was " << (EmptySetTraits<FromSetType>::isEmpty(m_fromSet) ? "" : " not ") << "null"
                  << " , toSet was " << (EmptySetTraits<ToSetType>::isEmpty(m_toSet) ? "" : " not ") << "null";
        }

        bValid = false;
      }
    }
    else
    {
      if(verboseOutput)
        errSstr << "\n\t* Neither set was null";

      // Check that the toSetIndices vector has the right size
      if( static_cast<SetPosition>(m_toSetIndicesVec.size()) != (stride() * m_fromSet->size()) )
      {
        if(verboseOutput)
        {
          errSstr << "\n\t* toSetIndices has the wrong size."
                  << "\n\t-- from set size is: " << m_fromSet->size()
                  << "\n\t-- constant stride is: " << stride()
                  << "\n\t-- expected relation size: " << (stride() * m_fromSet->size())
                  << "\n\t-- actual size: " << m_toSetIndicesVec.size()
          ;
        }
        bValid = false;
      }


      // Check that all elements of the toSetIndices vector point to valid set elements
      for(RelationVecConstIterator it = m_toSetIndicesVec.begin(), itEnd = m_toSetIndicesVec.end(); it != itEnd; ++it)
      {
        if( *it >= m_toSet->size() )
        {
          if(verboseOutput)
          {
            errSstr << "\n\t* toSetIndices had an out-of-range element."
                    << " -- value of element " << std::distance(m_toSetIndicesVec.begin(), it) << " was " << *it
                    << ". Max possible value should be " << m_toSet->size() << ".";
          }

          bValid = false;
        }
      }
    }


    if(verboseOutput)
    {
      std::stringstream sstr;

      if(bValid)
      {
        sstr << "(static,constant) Relation with stride " << stride() << " was valid." << std::endl;
      }
      else
      {
        sstr  << "Relation was NOT valid.\n"
              << errSstr.str()
              << std::endl;
      }

      sstr << "\n*** Detailed results of isValid on the relation.\n";
      if(m_fromSet)
        sstr << "\n** fromSet has size " << m_fromSet->size() << ": ";
      if(m_toSet)
        sstr << "\n** toSet has size " << m_toSet->size() << ": ";

      sstr << "\n** toSetIndices vec w/ size " << m_toSetIndicesVec.size() << ": ";
      std::copy(m_toSetIndicesVec.begin(), m_toSetIndicesVec.end(), std::ostream_iterator<SetPosition>(sstr, " "));

      std::cout << sstr.str() << std::endl;

    }

    return bValid;
  }



} // end namespace slam
} // end namespace asctoolkit

#endif // SLAM_STATIC_CONSTANT_RELATION_HPP_
