// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"
#include "axom/sidre.hpp"

using axom::sidre::DataStore;
using axom::sidre::Group;
using axom::sidre::View;
using axom::sidre::indexIsValid;
using axom::sidre::IndexType;
using axom::sidre::InvalidIndex;
using axom::sidre::InvalidName;
using axom::sidre::nameIsValid;

//-----Mock Serac Structs------------------------------------------------------

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

template <typename T, std::size_t Rows, std::size_t Cols>
class Tensor {
public:
    // Constructor to initialize all elements to a specific value
    Tensor(const T& value = T()) {
        for (auto& row : data) {
            row.fill(value);
        }
    }

    // Constructor to initialize from an initializer list
    Tensor(std::initializer_list<std::initializer_list<T>> init) {
        size_t row = 0;
        for (auto& rowData : init) {
            size_t col = 0;
            for (auto& element : rowData) {
                data[row][col++] = element;
            }
            row++;
        }
    }

    // Accessor for elements (read/write)
    T& at(size_t row, size_t col) {
        return data[row][col];
    }

    // Const accessor for elements (read-only)
    const T& at(size_t row, size_t col) const {
        return data[row][col];
    }

    // Get tensor dimensions
    std::pair<size_t, size_t> dimensions() const {
        return {Rows, Cols};
    }

    // Print the tensor for demonstration purposes
    void print() const {
        for (const auto& row : data) {
            for (const auto& element : row) {
                std::cout << element << " ";
            }
            std::cout << "\n";
        }
    }

private:
    std::array<std::array<T, Cols>, Rows> data;
};

//---------------------

template <>
bool check(const double& state, double value)
{
  SLIC_INFO(axom::fmt::format("{} == {}", state, value));
  return state == value;
}

template <>
void fill(double& state, double value)
{
  state = value;
}

//---------------------

struct StateOne {
  double x;
};

template <>
bool check(const StateOne& state, double value)
{
  SLIC_INFO(axom::fmt::format("{} == {}", state.x, value));
  return state.x == value;
}

template <>
void fill(StateOne& state, double value)
{
  state.x = value;
}

//---------------------

struct StateTwo {
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

struct StateThree {
  double x;
  double y;
  double z;
};

template <>
bool check(const StateThree& state, double value)
{
  return (state.x == value) && (state.y == (value + 1)) && (state.z == (value + 2));
}

template <>
void fill(StateThree& state, double value)
{
  state.x = value;
  state.y = value + 1;
  state.z = value + 2;
}

//---------------------

struct StateTensor {
  Tensor<double, 2, 2> t;
  double x;
};

template <>
bool check(const StateTensor& state, double value)
{
  return (state.t.at(0,0) == value) && (state.t.at(0,1) == (value + 1)) &&
         (state.t.at(1,0) == (value + 2)) && (state.t.at(1,1) == (value + 3)) &&
         (state.x == (value + 4));
}

template <>
void fill(StateTensor& state, double value)
{
  state.t = {{value, value + 1}, {value + 2, value + 3}};
  state.x = value + 4;
}

//---------------------

template <typename T>
using QuadratureData1D = axom::Array<T, 1>;
template <typename T>
using QuadratureData2D = axom::Array<T, 2>;

//------------------------------------------------------------------------------

template <typename T>
void test_user_defined_data()
{
  // populate data
  constexpr IndexType size = 10;
  QuadratureData1D<T> states(size, size);
  IndexType i = 0;
  for(auto& state: states){
    fill(state, i++);
  }

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

  // write data to datastore as bytes
  View* data_view = qd_group->createViewAndAllocate("states",
                                                    axom::sidre::UINT8_ID,
                                                    total_size);
  std::uint8_t* sidre_state_data = data_view->getData();
  memcpy(sidre_state_data, states.data(), static_cast<std::size_t>(total_size));

  // mess with data before overriding it back to original in sidre
  fill(states[0], 123);
  fill(states[size/2], 456);
  fill(states[size-1], 789);

  // Copy original data back over local data
  memcpy(states.data(), sidre_state_data, static_cast<std::size_t>(total_size));

  // Test data is back to original
  i = 0;
  for(auto& state: states){
    EXPECT_TRUE(check(state, i++));
  }
}

template <typename T>
void test_external_user_defined_data()
{
  // populate data
  constexpr IndexType size = 10;
  QuadratureData1D<T> states(size, size);
  IndexType i = 0;
  for(auto& state: states){
    fill(state, i++);
  }

  // Create datastore
  DataStore ds;
  Group* qd_group = ds.getRoot()->createGroup("quadraturedata");

  // get size
  auto num_states = static_cast<IndexType>(states.size());
  auto state_size = static_cast<IndexType>(sizeof(states[0]));
  auto total_size = num_states * state_size;
  auto num_uint8s = total_size / sizeof(axom::sidre::UINT8_ID);

  // write shape
  qd_group->createViewScalar("num_states", num_states);
  qd_group->createViewScalar("state_size", state_size);
  qd_group->createViewScalar("total_size", total_size);

  // Add states as an external buffer
  View* states_view = qd_group->createView("states");
  states_view->setExternalDataPtr(axom::sidre::UINT8_ID, num_uint8s, states.data());

  // Save the array data in to a file
  std::string filename = "sidre_external_quadraturedata";
  qd_group->save(filename);

  // Create new array to fill with saved data
  QuadratureData1D<T> saved_states(size, size);
  for(auto& state: saved_states){
    fill(state, -1);
  }

  // Load data into new array
  Group* new_qd_group = ds.getRoot()->createGroup("new_quadraturedata");
  new_qd_group->load(filename);
  EXPECT_TRUE(new_qd_group->hasView("states"));
  View* new_states_view = new_qd_group->getView("states");
  new_states_view->setExternalDataPtr(saved_states.data());
  new_qd_group->loadExternalData(filename);

  // Test data is back to original
  i = 0;
  for(auto& state: saved_states){
    EXPECT_TRUE(check(state, i++));
  }
}

//------------------------------------------------------------------------------

TEST(sidre, QD_double_readandwrite)
{
  test_user_defined_data<double>();
}

TEST(sidre, QD_StateOne_readandwrite)
{
  test_user_defined_data<StateOne>();
}

TEST(sidre, QD_StateTwo_readandwrite)
{
  test_user_defined_data<StateTwo>();
}

TEST(sidre, QD_StateThree_readandwrite)
{
  test_user_defined_data<StateThree>();
}

TEST(sidre, QD_StateTensor_readandwrite)
{
  test_user_defined_data<StateTensor>();
}

//------------------------------------------------------------------------------

TEST(sidre, QD_double_external_readandwrite)
{
  test_external_user_defined_data<double>();
}

TEST(sidre, QD_StateOne_external_readandwrite)
{
  test_external_user_defined_data<StateOne>();
}

//----------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}

