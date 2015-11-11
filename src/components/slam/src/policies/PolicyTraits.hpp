#ifndef SLAM_POLICY_TRAITS_H_
#define SLAM_POLICY_TRAITS_H_


#include "slam/SizePolicies.hpp"
#include "slam/StridePolicies.hpp"

namespace asctoolkit {
namespace slam {
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
} // end namespace slam
} // end namespace asctoolkit

#endif // SLAM_POLICY_TRAITS_H_
