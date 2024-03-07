// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file core_array_perf.cpp
 *  \brief This example illustrates performance of array classes.
 */

// Axom includes
#include "axom/core/StackArray.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/ArrayIndexer.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/runtime_policy.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/core/utilities/Timer.hpp"
#include "axom/core/utilities/AnnotationMacros.hpp"

#include "axom/fmt/format.h"
#include "axom/CLI11.hpp"

// C/C++ includes
#include <iostream>
#include <vector>

///////////////////////////////////////////////////////////////
/// Struct to parse and store the input parameters
struct Input
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;

  // Array shape.
  std::vector<axom::IndexType> shape;
  // Array stride order (same length as shape)
  std::vector<std::uint32_t> slowestDirections;
  axom::ArrayStrideOrder strideOrder = axom::ArrayStrideOrder::ARBITRARY;

  RuntimePolicy policy = RuntimePolicy::seq;

  bool checkResults = false;

private:
  bool _verboseOutput {false};
  const std::map<std::string, int> strideValidator {
    {"row", int(axom::ArrayStrideOrder::ROW)},
    {"col", int(axom::ArrayStrideOrder::COLUMN)}};

public:
  bool isVerbose() const { return _verboseOutput; }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-p, --policy", policy)
      ->description("Set runtime policy for test")
      ->capture_default_str()
      ->transform(
        axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

    app.add_flag("-v,--verbose,!--no-verbose", _verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_flag("-c,--check-results,!--no-check-results", checkResults)
      ->description("Enable/disable checking against expected results")
      ->capture_default_str();

    app.add_option("--shape", shape)->description("Array shape")->expected(1, 4);

    auto* strideOrderOption =
      app.add_option("--strideOrder", strideOrder)
        ->description("Stride order (must be as long as sizes.")
        ->transform(axom::CLI::CheckedTransformer(strideValidator))
        ->expected(1, 4);

    app.add_option("--slowestDirections", slowestDirections)
      ->description(
        "Stride directions, from slowest to fastest.  Must be same length as "
        "shape.")
      ->excludes(strideOrderOption);

    app.get_formatter()->column_width(60);

    app.parse(argc, argv);

    if (shape.empty())
    {
      std::cerr << "You must specify shape (1-4 integers)." << std::endl;
      std::abort();
    }

    /*
      If slowestDirections is specified, it must match length of shape.
      If neither is specified, default to row-major ordering.
    */
    if(!slowestDirections.empty())
    {
      if (slowestDirections.size() != shape.size())
      {
        std::cerr << "slowestDimension size (" << slowestDirections.size()
                  << ") must match shape size (" << shape.size() << ")."
                  << std::endl;
        std::abort();
      }
    }
    if (slowestDirections.empty() &&
        strideOrder == axom::ArrayStrideOrder::ARBITRARY)
    {
      strideOrder = axom::ArrayStrideOrder::ROW;
    }
  }
};

Input params;

template <typename T>
std::string array_to_string(const T* dataPtr, int count)
{
  std::ostringstream os;
  os << '[';
  for(int i = 0; i < count; ++i)
  {
    os << dataPtr[i] << (i < count - 1 ? ',' : ']');
  }
  return os.str();
}

/*!
  @brief Time the sequential access of every element of an array.

  This is the fastest we expect to visit every element.
*/
template <int DIM>
void runTest_sequentialAccess(axom::ArrayView<double, DIM>& array)
{
  auto count = array.size();
  auto* ptr = array.data();
  for(axom::IndexType i = 0; i < count; ++i)
  {
    ptr[i] += 1;
  }
}

/*!
  @brief Time the row-major access of every element of an array.
*/
template <int DIM>
void runTest_rowMajorAccess(axom::ArrayView<double, DIM>& array);

/*!
  @brief Time the column-major access of every element of an array.
*/
template <int DIM>
void runTest_rowMajorAccess(axom::ArrayView<double, DIM>& array);

/*!
  @brief Time the access of every elemenent in array in
  a dynamic loop ordering.

  The dynamic order should match the sequential order,
  Any performance difference is due to overhead of dynamic
  nesting of the loops.
*/
template <int DIM>
void runTest_dynamicAccess(axom::ArrayView<double, DIM>& array);

void runTest_rowMajorAccess(axom::ArrayView<double, 1>& array)
{
  AXOM_PERF_MARK_FUNCTION("rowMajorAccess-1D");
  const auto& shape = array.shape();
  for(axom::IndexType i = 0; i < shape[0]; ++i)
  {
    array[i] += 10;
  }
}

void runTest_rowMajorAccess(axom::ArrayView<double, 2>& array)
{
  AXOM_PERF_MARK_FUNCTION("rowMajorAccess-2D");
  const auto& shape = array.shape();
  for(axom::IndexType i = 0; i < shape[0]; ++i)
  {
    for(axom::IndexType j = 0; j < shape[1]; ++j)
    {
      array(i, j) += 10;
    }
  }
}

void runTest_rowMajorAccess(axom::ArrayView<double, 3>& array)
{
  AXOM_PERF_MARK_FUNCTION("rowMajorAccess-3D");
  const auto& shape = array.shape();
  for(axom::IndexType i = 0; i < shape[0]; ++i)
  {
    for(axom::IndexType j = 0; j < shape[1]; ++j)
    {
      for(axom::IndexType k = 0; k < shape[2]; ++k)
      {
        array(i, j, k) += 10;
      }
    }
  }
}

void runTest_rowMajorAccess(axom::ArrayView<double, 4>& array)
{
  AXOM_PERF_MARK_FUNCTION("rowMajorAccess-4D");
  const auto& shape = array.shape();
  for(axom::IndexType i = 0; i < shape[0]; ++i)
  {
    for(axom::IndexType j = 0; j < shape[1]; ++j)
    {
      for(axom::IndexType k = 0; k < shape[2]; ++k)
      {
        for(axom::IndexType l = 0; l < shape[3]; ++l)
        {
          array(i, j, k, l) += 10;
        }
      }
    }
  }
}

void runTest_columnMajorAccess(axom::ArrayView<double, 1>& array)
{
  AXOM_PERF_MARK_FUNCTION("columnMajorAccess-1D");
  const auto& shape = array.shape();
  for(axom::IndexType i = 0; i < shape[0]; ++i)
  {
    array[i] += 100;
  }
}

void runTest_columnMajorAccess(axom::ArrayView<double, 2>& array)
{
  AXOM_PERF_MARK_FUNCTION("columnMajorAccess-2D");
  const auto& shape = array.shape();
  for(axom::IndexType j = 0; j < shape[1]; ++j)
  {
    for(axom::IndexType i = 0; i < shape[0]; ++i)
    {
      array(i, j) += 100;
    }
  }
}

void runTest_columnMajorAccess(axom::ArrayView<double, 3>& array)
{
  AXOM_PERF_MARK_FUNCTION("columnMajorAccess-3D");
  const auto& shape = array.shape();
  for(axom::IndexType k = 0; k < shape[2]; ++k)
  {
    for(axom::IndexType j = 0; j < shape[1]; ++j)
    {
      for(axom::IndexType i = 0; i < shape[0]; ++i)
      {
        array(i, j, k) += 100;
      }
    }
  }
}

void runTest_columnMajorAccess(axom::ArrayView<double, 4>& array)
{
  AXOM_PERF_MARK_FUNCTION("columnMajorAccess-4D");
  const auto& shape = array.shape();
  for(axom::IndexType l = 0; l < shape[3]; ++l)
  {
    for(axom::IndexType k = 0; k < shape[2]; ++k)
    {
      for(axom::IndexType j = 0; j < shape[1]; ++j)
      {
        for(axom::IndexType i = 0; i < shape[0]; ++i)
        {
          array(i, j, k, l) += 100;
        }
      }
    }
  }
}

void runTest_dynamicAccess(axom::ArrayView<double, 1>& array)
{
  AXOM_PERF_MARK_FUNCTION("dynamicAccess-1D");
  const auto& shape = array.shape();
  for(axom::IndexType i = 0; i < shape[0]; ++i)
  {
    array[i] += 1000;
  }
}

void runTest_dynamicAccess(axom::ArrayView<double, 2>& array)
{
  AXOM_PERF_MARK_FUNCTION("dynamicAccess-2D");
  const auto& shape = array.shape();
  const auto& indexer = array.indexer();
  const auto& slowestDirs = indexer.slowestDirs();
  axom::StackArray<axom::IndexType, 2> idx;
  axom::IndexType& m = idx[slowestDirs[0]];
  axom::IndexType& n = idx[slowestDirs[1]];
  const axom::StackArray<axom::IndexType, 2> ends {shape[slowestDirs[0]],
                                                   shape[slowestDirs[1]]};
  for(m = 0; m < ends[0]; ++m)
  {
    for(n = 0; n < ends[1]; ++n)
    {
      array[idx] += 1000;
    }
  }
}

void runTest_dynamicAccess(axom::ArrayView<double, 3>& array)
{
  AXOM_PERF_MARK_FUNCTION("dynamicAccess-3D");
  const auto& shape = array.shape();
  const auto& indexer = array.indexer();
  const auto& slowestDirs = indexer.slowestDirs();
  axom::StackArray<axom::IndexType, 3> idx;
  axom::IndexType& m = idx[slowestDirs[0]];
  axom::IndexType& n = idx[slowestDirs[1]];
  axom::IndexType& o = idx[slowestDirs[2]];
  const axom::StackArray<axom::IndexType, 3> ends {shape[slowestDirs[0]],
                                                   shape[slowestDirs[1]],
                                                   shape[slowestDirs[2]]};
  for(m = 0; m < ends[0]; ++m)
  {
    for(n = 0; n < ends[1]; ++n)
    {
      for(o = 0; o < ends[2]; ++o)
      {
        array[idx] += 1000;
      }
    }
  }
}

void runTest_dynamicAccess(axom::ArrayView<double, 4>& array)
{
  AXOM_PERF_MARK_FUNCTION("dynamicAccess-4D");
  const auto& shape = array.shape();
  const auto& indexer = array.indexer();
  const auto& slowestDirs = indexer.slowestDirs();
  axom::StackArray<axom::IndexType, 4> idx;
  axom::IndexType& m = idx[slowestDirs[0]];
  axom::IndexType& n = idx[slowestDirs[1]];
  axom::IndexType& o = idx[slowestDirs[2]];
  axom::IndexType& p = idx[slowestDirs[3]];
  const axom::StackArray<axom::IndexType, 4> ends {shape[slowestDirs[0]],
                                                   shape[slowestDirs[1]],
                                                   shape[slowestDirs[2]],
                                                   shape[slowestDirs[3]]};
  for(m = 0; m < ends[0]; ++m)
  {
    for(n = 0; n < ends[1]; ++n)
    {
      for(o = 0; o < ends[2]; ++o)
      {
        for(p = 0; p < ends[3]; ++p)
        {
          array[idx] += 1000;
        }
      }
    }
  }
}

/*!
  @brief Return an array for testing, dimension DIM,
  sized and ordered according to params values.
*/
template <int DIM>
axom::Array<double, DIM> makeArray()
{
  assert(DIM <= params.shape.size());

  axom::StackArray<axom::IndexType, DIM> shape;
  for(int d = 0; d < DIM; ++d)
  {
    shape[d] = params.shape[d];
  }

  // Until indexer isused, array will have row-major order.
  assert(params.slowestDirections.empty());
  assert(int(params.strideOrder) & int(axom::ArrayStrideOrder::ROW));

  axom::ArrayIndexer<axom::IndexType, DIM> indexer;
  if(!params.slowestDirections.empty())
  {
    axom::StackArray<std::uint16_t, DIM> slowestDirections;
    for(int d = 0; d < DIM; ++d)
    {
      slowestDirections[d] = params.slowestDirections[d];
    }
    indexer.initializeShape(shape, slowestDirections);
  }
  else
  {
    indexer.initializeShape(shape, params.strideOrder);
  }

  return axom::Array<double, DIM>(shape);
}

/*!
  @brief Return an array for testing, dimension DIM,
  sized and ordered according to params values.

  This method allocates a 1D array and puts a muldimensional
  view on the data.  The view supports arbitrary ordering.
*/
template <int DIM>
void makeArray(axom::Array<double> & ar, axom::ArrayView<double, DIM> & view)
{
  assert(DIM <= params.shape.size());

  axom::StackArray<axom::IndexType, DIM> shape;
  for(int d = 0; d < DIM; ++d)
  {
    shape[d] = params.shape[d];
  }

  axom::IndexType product = shape[0];
  for(int d = 1; d < DIM; ++d)
  {
    product *= shape[d];
  }

  ar.resize(product);

  axom::ArrayIndexer<axom::IndexType, DIM> indexer;
  if(!params.slowestDirections.empty())
  {
    axom::StackArray<std::uint16_t, DIM> slowestDirections;
    for(int d = 0; d < DIM; ++d)
    {
      slowestDirections[d] = params.slowestDirections[d];
    }
    indexer.initializeShape(shape, slowestDirections);
  }
  else
  {
    indexer.initializeShape(shape, params.strideOrder);
  }
  view = axom::ArrayView<double, DIM>(ar.data(), shape, indexer.strides());

  return;
}

/*!
  @brief Run test with array shape of the first DIM values in
  params.shape.
*/
template <int DIM>
void runTest_dim()
{
  // Use ArrayView to test, because Array doesn't support
  // arbitrary ordering (yet).
#if 0
  axom::Array<double, DIM> arrayMd = makeArray<DIM>();
  axom::ArrayView<double, DIM> array = arrayMd.view();
#else
  axom::Array<double> array1D;
  axom::ArrayView<double, DIM> array;
  makeArray<DIM>(array1D, array);
#endif

  auto count = array.size();
  auto* ptr = array.data();
  for(axom::IndexType i = 0; i < count; ++i)
  {
    ptr[i] = double(i*1000000);
  }

  axom::utilities::Timer sequentialTimer(false);
  sequentialTimer.start();
  runTest_sequentialAccess(array);
  sequentialTimer.stop();
  std::cout << "Sequential time   " << sequentialTimer.elapsedTimeInSec() << " seconds, base" << std::endl;

  axom::utilities::Timer rowMajorTimer(false);
  rowMajorTimer.start();
  runTest_rowMajorAccess(array);
  rowMajorTimer.stop();
  std::cout << "Row-major time    " << rowMajorTimer.elapsedTimeInSec() << " seconds, "
            << std::setprecision(3) << rowMajorTimer.elapsedTimeInSec()/sequentialTimer.elapsedTimeInSec() << 'x' << std::endl;

  axom::utilities::Timer columnMajorTimer(false);
  columnMajorTimer.start();
  runTest_columnMajorAccess(array);
  columnMajorTimer.stop();
  std::cout << "Column-major time " << columnMajorTimer.elapsedTimeInSec() << " seconds, "
            << std::setprecision(3) << columnMajorTimer.elapsedTimeInSec()/sequentialTimer.elapsedTimeInSec() << 'x'  << std::endl;

  axom::utilities::Timer dynamicTimer(false);
  dynamicTimer.start();
  runTest_dynamicAccess(array);
  dynamicTimer.stop();
  std::cout << "Dynamic time      " << dynamicTimer.elapsedTimeInSec() << " seconds, "
            << std::setprecision(3) << dynamicTimer.elapsedTimeInSec()/sequentialTimer.elapsedTimeInSec() << 'x'  << std::endl;

}

/*!
  @brief Run test based on dimension specified in params.
*/
void runTest()
{
  switch(params.shape.size())
  {
  case 1:
    runTest_dim<1>();
    break;
  case 2:
    runTest_dim<2>();
    break;
  case 3:
    runTest_dim<3>();
    break;
  case 4:
    runTest_dim<4>();
    break;
  }
}

int main(int argc, char** argv)
{
  axom::CLI::App app {"Driver for array indexing performace tests"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    auto retval = app.exit(e);
    exit(retval);
  }

  runTest();

  return 0;
}
