// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/sidre.hpp"

using axom::sidre::DataStore;
using axom::sidre::Group;
using axom::sidre::IndexType;
using axom::sidre::View;

//-----Mock Serac Structs------------------------------------------------------

/**
 * @brief Arbitrary-rank tensor class
 * @tparam T The type stored at each index
 * @tparam n The dimensions of the tensor
 */
template <typename T, int... n>
struct tensor;

template <typename T, int m, int... n>
struct tensor<T, m, n...>
{
  template <typename i_type>
  constexpr auto& operator()(i_type i)
  {
    return data[i];
  }
  template <typename i_type>
  constexpr auto& operator()(i_type i) const
  {
    return data[i];
  }
  template <typename i_type, typename... jklm_type>
  constexpr auto& operator()(i_type i, jklm_type... jklm)
  {
    return data[i](jklm...);
  }
  template <typename i_type, typename... jklm_type>
  constexpr auto& operator()(i_type i, jklm_type... jklm) const
  {
    return data[i](jklm...);
  }

  constexpr auto& operator[](int i) { return data[i]; }
  constexpr const auto& operator[](int i) const { return data[i]; }

  tensor<T, n...> data[m];
};

template <typename T, int m>
struct tensor<T, m>
{
  template <typename i_type>
  constexpr auto& operator()(i_type i)
  {
    return data[i];
  }
  template <typename i_type>
  constexpr auto& operator()(i_type i) const
  {
    return data[i];
  }
  constexpr auto& operator[](int i) { return data[i]; }
  constexpr const auto& operator[](int i) const { return data[i]; }

  template <int last_dimension = m,
            typename = typename std::enable_if<last_dimension == 1>::type>
  constexpr operator T()
  {
    return data[0];
  }

  template <int last_dimension = m,
            typename = typename std::enable_if<last_dimension == 1>::type>
  constexpr operator T() const
  {
    return data[0];
  }

  T data[m];
};

//---------------------

template <typename T>
bool check(const T&, double)
{
  SLIC_ERROR("You didn't implement a specialization for check");
  return false;
}

template <typename T>
void fill(T&, double)
{
  SLIC_ERROR("You didn't implement a specialization for fill");
}

//---------------------

template <>
bool check(const double& state, double value)
{
  return state == value;
}

template <>
void fill(double& state, double value)
{
  state = value;
}

//---------------------

struct StateOne
{
  double x;
};

template <>
bool check(const StateOne& state, double value)
{
  return state.x == value;
}

template <>
void fill(StateOne& state, double value)
{
  state.x = value;
}

//---------------------

struct StateTwo
{
  double x;
  double y;
};

template <>
bool check(const StateTwo& state, double value)
{
  return (state.x == value) && (state.y == (value + 1));
}

template <>
void fill(StateTwo& state, double value)
{
  state.x = value;
  state.y = value + 1;
}

//---------------------

struct StateThree
{
  double x;
  double y;
  double z;
};

template <>
bool check(const StateThree& state, double value)
{
  return (state.x == value) && (state.y == (value + 1)) &&
    (state.z == (value + 2));
}

template <>
void fill(StateThree& state, double value)
{
  state.x = value;
  state.y = value + 1;
  state.z = value + 2;
}

//---------------------

struct StateTensor
{
  //tensor<double, 2, 2> t;
  tensor<double, 2> t;
  double x;
};

template <>
bool check(const StateTensor& state, double value)
{
  // SLIC_INFO(axom::fmt::format("check({}, {}, {}, {}), {}",
  //                             state.t(0, 0),
  //                             state.t(0, 1),
  //                             state.t(1, 0),
  //                             state.t(1, 1),
  //                             state.x));

  // return (state.t(0, 0) == value) && (state.t(0, 1) == (value + 1)) &&
  //   (state.t(1, 0) == (value + 2)) && (state.t(1, 1) == (value + 3)) &&
  //   (state.x == (value + 4));

  SLIC_INFO(
    axom::fmt::format("check({}, {}), {}", state.t(0), state.t(1), state.x));

  return (state.t(0) == value) && (state.t(1) == (value + 1)) &&
    (state.x == (value + 4));
}

template <>
void fill(StateTensor& state, double value)
{
  // state.t(0,0) = value;
  // state.t(0,1) = value + 1;
  // state.t(1,0) = value + 2;
  // state.t(1,1) = value + 3;
  state.t(0) = value;
  state.t(1) = value + 1;
  state.x = value + 4;
}

//---------------------

template <typename T, int DIM>
void fill_array(T, int, bool)
{
  SLIC_ERROR("You didn't implement the fill array of this type.");
}

template <typename T, int DIM>
void check_array(T)
{
  SLIC_ERROR("You didn't implement the check array of this type.");
}

template <typename T, int DIM>
void fill_array(axom::Array<T, DIM>& states, double fill_value, bool increment)
{
  int count = 0;
  for(auto& state : states)
  {
    fill(state, fill_value);
    if(increment)
    {
      fill_value++;
    }
    count++;
  }
  SLIC_INFO(
    axom::fmt::format("filled {} states with end value {}", count, fill_value));
}

template <typename T, int DIM>
void check_array(axom::Array<T, DIM> states)
{
  int count = 0;
  double check_value = 0;
  for(auto& state : states)
  {
    EXPECT_TRUE(check(state, check_value++));
    count++;
  }
  SLIC_INFO(
    axom::fmt::format("checked {} states with end value {}", count, check_value));
}

//------------------------------------------------------------------------------

template <typename T, int DIM>
void test_user_defined_data()
{
  // populate data
  constexpr IndexType size = 10;
  axom::Array<T, DIM> states(size, size);
  fill_array(states, 0, true);

  // Create datastore
  DataStore ds;
  Group* qd_group = ds.getRoot()->createGroup("quadraturedata");

  // get size
  auto num_states = static_cast<IndexType>(states.size());
  auto state_size = static_cast<IndexType>(sizeof(states[0]));
  auto total_size = num_states * state_size;

  // write shape
  qd_group->createViewScalar("num_states", num_states);
  qd_group->createViewScalar("state_size", state_size);
  qd_group->createViewScalar("total_size", total_size);
  SLIC_INFO(axom::fmt::format("Num of States={}", num_states));
  SLIC_INFO(axom::fmt::format("State Size={}", state_size));
  SLIC_INFO(axom::fmt::format("Total Size={}", total_size));
  SLIC_INFO(axom::fmt::format("Double Size={}", sizeof(double)));
  SLIC_INFO(axom::fmt::format("Total Size/Double Size={}",
                              total_size / sizeof(double)));

  // write data to datastore as bytes
  View* data_view =
    qd_group->createViewAndAllocate("states", axom::sidre::UINT8_ID, total_size);
  std::uint8_t* sidre_state_data = data_view->getData();
  memcpy(sidre_state_data, states.data(), static_cast<std::size_t>(total_size));

  // mess with data before overriding it back to original in sidre
  fill_array(states, -1, false);

  // Copy original data back over local data
  memcpy(states.data(), sidre_state_data, static_cast<std::size_t>(total_size));

  // Test data is back to original
  check_array(states);
}

template <typename T, int DIM>
void test_external_user_defined_data()
{
  // populate data
  constexpr IndexType size = 10;
  axom::Array<T, DIM> states(size, size);
  fill_array(states, 0, true);

  // Create datastore
  DataStore ds;
  Group* qd_group = ds.getRoot()->createGroup("quadraturedata");

  // get size
  auto num_states = static_cast<IndexType>(states.size());
  auto state_size = static_cast<IndexType>(sizeof(states[0]));
  auto total_size = num_states * state_size;
  auto num_uint8s = total_size / sizeof(std::uint8_t);

  // write shape
  qd_group->createViewScalar("num_states", num_states);
  qd_group->createViewScalar("state_size", state_size);
  qd_group->createViewScalar("total_size", total_size);

  // Add states as an external buffer
  View* states_view = qd_group->createView("states");
  states_view->setExternalDataPtr(axom::sidre::UINT8_ID,
                                  num_uint8s,
                                  states.data());

  // Save the array data in to a file
  std::string filename = "sidre_external_quadraturedata";
  qd_group->save(filename);

  // Create new array to fill with saved data
  axom::Array<T, DIM> saved_states(size, size);
  fill_array(saved_states, -1, false);

  // Load data into new array
  Group* new_qd_group = ds.getRoot()->createGroup("new_quadraturedata");
  new_qd_group->load(filename);
  EXPECT_TRUE(new_qd_group->hasView("states"));
  View* new_states_view = new_qd_group->getView("states");
  new_states_view->setExternalDataPtr(saved_states.data());
  new_qd_group->loadExternalData(filename);

  // Test data was read into the new location and overwrote the bad data
  check_array(saved_states);
}

//------------------------------------------------------------------------------

TEST(sidre, OneD_double_readandwrite) { test_user_defined_data<double, 1>(); }

// TEST(sidre, OneD_StateOne_readandwrite) { test_user_defined_data<StateOne, 1>(); }

// TEST(sidre, OneD_StateTwo_readandwrite) { test_user_defined_data<StateTwo, 1>(); }

// TEST(sidre, OneD_StateThree_readandwrite) { test_user_defined_data<StateThree, 1>(); }

// TEST(sidre, OneD_StateTensor_readandwrite) { test_user_defined_data<StateTensor, 1>(); }

// //-------------------------

// TEST(sidre, TwoD_double_readandwrite) { test_user_defined_data<double, 2>(); }

// TEST(sidre, TwoD_StateOne_readandwrite) { test_user_defined_data<StateOne, 2>(); }

// TEST(sidre, TwoD_StateTwo_readandwrite) { test_user_defined_data<StateTwo, 2>(); }

// TEST(sidre, TwoD_StateThree_readandwrite) { test_user_defined_data<StateThree, 2>(); }

TEST(sidre, TwoD_StateTensor_readandwrite)
{
  test_user_defined_data<StateTensor, 2>();
}

//------------------------------------------------------------------------------

// TEST(sidre, OneD_double_external_readandwrite) { test_external_user_defined_data<double, 1>(); }

// TEST(sidre, OneD_StateOne_external_readandwrite) { test_external_user_defined_data<StateOne, 1>(); }

// TEST(sidre, OneD_StateTwo_external_readandwrite) { test_external_user_defined_data<StateTwo, 1>(); }

// TEST(sidre, OneD_StateThree_external_readandwrite) { test_external_user_defined_data<StateThree, 1>(); }

// TEST(sidre, OneD_StateTensor_external_readandwrite) { test_external_user_defined_data<StateTensor, 1>(); }

//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
