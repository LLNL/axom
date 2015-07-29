/**
 * \file IndirectionSet.hpp
 *
 * \brief Basic API for a set of entities in a simulation
 */

#ifndef MESHAPI_INDIRECTION_SET_H_
#define MESHAPI_INDIRECTION_SET_H_

#include <cstddef>
#include <vector>

#include "meshapi/OrderedSet.hpp"

namespace asctoolkit {
namespace meshapi {


/**
 * \class IndirectionSet
 *
 * \brief An indexed set (a tuple) of entities in a simulation
 *
 * A container class for a set of entities in a simulation. Each entity has an index.
 *
 * Below is an initial implementation for a set with explicit indexes (encoded here using a vector).
 */
  class ArrayIndirectionSet : public OrderedSet< policies::RuntimeSizeHolder<Set::PositionType>
                                               , policies::ZeroOffset<Set::PositionType>
                                               , policies::StrideOne<Set::PositionType>
                                               , policies::ArrayIndirection<Set::PositionType, Set::ElementType>
                                               >
  {
  private:
      typedef OrderedSet< policies::RuntimeSizeHolder<Set::PositionType>
      , policies::ZeroOffset<Set::PositionType>
      , policies::StrideOne<Set::PositionType>
      , policies::ArrayIndirection<Set::PositionType, Set::ElementType>
      >  OrderedSetType;

  public:
      typedef typename OrderedSetType::PositionType PositionType;
      typedef typename OrderedSetType::IndexType IndexType;
      typedef typename OrderedSetType::ElementType ElementType;

  private:
      static const PositionType DEFAULT_SIZE = OrderedSetType::SizePolicyType::DEFAULT_VALUE;
      static const PositionType DEFAULT_OFFSET = OrderedSetType::OffsetPolicyType::DEFAULT_VALUE;
      static const PositionType DEFAULT_STRIDE = OrderedSetType::StridePolicyType::DEFAULT_VALUE;

  public:
      ArrayIndirectionSet (PositionType size = DEFAULT_SIZE)
          : OrderedSetType(size, DEFAULT_OFFSET, DEFAULT_STRIDE) {}
    ~ArrayIndirectionSet () {}
  };

  class VectorIndirectionSet : public OrderedSet< policies::RuntimeSizeHolder<Set::PositionType>
                                               , policies::ZeroOffset<Set::PositionType>
                                               , policies::StrideOne<Set::PositionType>
                                               , policies::STLVectorIndirection<Set::PositionType, Set::ElementType>
                                               // add parent subset ?
                                               >
  {
  private:
      typedef OrderedSet< policies::RuntimeSizeHolder<Set::PositionType>
      , policies::ZeroOffset<Set::PositionType>
      , policies::StrideOne<Set::PositionType>
      , policies::STLVectorIndirection<Set::PositionType, Set::ElementType>
      >  OrderedSetType;

      typedef OrderedSet::IndirectionPolicyType IndirectionPolicyType;

  public:
      typedef typename OrderedSetType::PositionType PositionType;
      typedef typename OrderedSetType::IndexType IndexType;
      typedef typename OrderedSetType::ElementType ElementType;

      typedef typename IndirectionPolicyType::VectorType ArrType;
  private:
      static const PositionType DEFAULT_SIZE = OrderedSetType::SizePolicyType::DEFAULT_VALUE;
      static const PositionType DEFAULT_OFFSET = OrderedSetType::OffsetPolicyType::DEFAULT_VALUE;
      static const PositionType DEFAULT_STRIDE = OrderedSetType::StridePolicyType::DEFAULT_VALUE;

  public:
      VectorIndirectionSet (PositionType size = DEFAULT_SIZE)
          : OrderedSetType(size, DEFAULT_OFFSET, DEFAULT_STRIDE) {}

    ~VectorIndirectionSet () {}
  };


} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_INDIRECTION_SET_H_
