#ifndef MESHAPI_POLICY_TRAITS_H_
#define MESHAPI_POLICY_TRAITS_H_


#include "meshapi/SizePolicies.hpp"
#include "meshapi/StridePolicies.hpp"

namespace asctoolkit {
namespace meshapi {
namespace policies {

  template<typename StridePolicyType, typename IntType, IntType VAL> struct StrideToSize;

  template<> struct StrideToSize< RuntimeStrideHolder< Set::PositionType >
                                , Set::PositionType
                                , RuntimeStrideHolder< Set::PositionType >::DEFAULT_VALUE >
  {
      typedef RuntimeSizeHolder<typename Set::PositionType> SizeType;
  };

  template<Set::PositionType VAL> struct StrideToSize< CompileTimeStrideHolder<Set::PositionType, VAL>, Set::PositionType, VAL >
  {
      typedef CompileTimeSizeHolder<Set::PositionType, VAL> SizeType;
  };

  template<> struct StrideToSize< StrideOne<Set::PositionType>, Set::PositionType,  StrideOne<Set::PositionType>::DEFAULT_VALUE >
  {
      typedef CompileTimeSizeHolder<Set::PositionType, StrideOne<Set::PositionType  >::DEFAULT_VALUE > SizeType;
  };


} // end namespace policies
} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_POLICY_TRAITS_H_
