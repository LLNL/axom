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
#include "axom/core/execution/nested_for_exec.hpp"
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
struct InputParams
{
public:
  using RuntimePolicy = axom::runtime_policy::Policy;

  // Array shape.
  std::vector<axom::IndexType> shape;
  axom::IndexType ghostWidth = 1;
  std::vector<axom::IndexType> paddedShape;
  std::vector<axom::IndexType> idxBegin;
  std::vector<axom::IndexType> idxEnd;
  axom::IndexType realSize;
  axom::IndexType paddedSize;

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

    app.add_option("--ghost", ghostWidth)->description("Ghost width");

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

    if(shape.empty())
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
      if(slowestDirections.size() != shape.size())
      {
        std::cerr << "slowestDimension size (" << slowestDirections.size()
                  << ") must match shape size (" << shape.size() << ")."
                  << std::endl;
        std::abort();
      }
    }
    if(slowestDirections.empty() &&
       strideOrder == axom::ArrayStrideOrder::ARBITRARY)
    {
      strideOrder = axom::ArrayStrideOrder::ROW;
    }

    //
    // Dependent data
    //

    paddedShape.resize(shape.size());
    idxBegin.resize(shape.size());
    idxEnd.resize(shape.size());
    for(size_t i = 0; i < shape.size(); ++i)
    {
      paddedShape[i] = shape[i] + 2 * ghostWidth;
      idxBegin[i] = ghostWidth;
      idxEnd[i] = idxBegin[i] + shape[i];
    }

    realSize = shape[0];
    paddedSize = paddedShape[0];
    for(size_t i = 1; i < shape.size(); ++i)
    {
      realSize *= shape[i];
      paddedSize *= paddedShape[i];
    }
  }
};

InputParams params;

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
#ifdef AXOM_USE_RAJA
  for(axom::IndexType i = 0; i < count; ++i)
  {
    ptr[i] += 1;
  }
#else
  for(axom::IndexType i = 0; i < count; ++i)
  {
    ptr[i] += 1;
  }
#endif
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

//
// Row-major access tests
//

void runTest_rowMajorAccess(axom::ArrayView<double, 1>& array)
{
  AXOM_PERF_MARK_FUNCTION("rowMajorAccess-1D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    array[i] += 10;
  }
#else
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    array[i] += 10;
  }
#endif
}

void runTest_rowMajorAccess(axom::ArrayView<double, 2>& array)
{
  AXOM_PERF_MARK_FUNCTION("rowMajorAccess-2D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
    {
      array(i, j) += 10;
    }
  }
#else
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
    {
      array(i, j) += 10;
    }
  }
#endif
}

void runTest_rowMajorAccess(axom::ArrayView<double, 3>& array)
{
  AXOM_PERF_MARK_FUNCTION("rowMajorAccess-3D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
    {
      for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
      {
        array(i, j, k) += 10;
      }
    }
  }
#else
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
    {
      for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
      {
        array(i, j, k) += 10;
      }
    }
  }
#endif
}

void runTest_rowMajorAccess(axom::ArrayView<double, 4>& array)
{
  AXOM_PERF_MARK_FUNCTION("rowMajorAccess-4D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
    {
      for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
      {
        for(axom::IndexType l = idxBegin[3]; l < idxEnd[3]; ++l)
        {
          array(i, j, k, l) += 10;
        }
      }
    }
  }
#else
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
    {
      for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
      {
        for(axom::IndexType l = idxBegin[3]; l < idxEnd[3]; ++l)
        {
          array(i, j, k, l) += 10;
        }
      }
    }
  }
#endif
}

//
// Colunn-major access tests
//

void runTest_columnMajorAccess(axom::ArrayView<double, 1>& array)
{
  AXOM_PERF_MARK_FUNCTION("columnMajorAccess-1D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    array[i] += 100;
  }
#else
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    array[i] += 100;
  }
#endif
}

void runTest_columnMajorAccess(axom::ArrayView<double, 2>& array)
{
  AXOM_PERF_MARK_FUNCTION("columnMajorAccess-2D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
  for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
  {
    for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
    {
      array(i, j) += 100;
    }
  }
#else
  for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
  {
    for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
    {
      array(i, j) += 100;
    }
  }
#endif
}

void runTest_columnMajorAccess(axom::ArrayView<double, 3>& array)
{
  AXOM_PERF_MARK_FUNCTION("columnMajorAccess-3D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
  for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
  {
    for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
    {
      for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
      {
        array(i, j, k) += 100;
      }
    }
  }
#else
  for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
  {
    for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
    {
      for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
      {
        array(i, j, k) += 100;
      }
    }
  }
#endif
}

void runTest_columnMajorAccess(axom::ArrayView<double, 4>& array)
{
  AXOM_PERF_MARK_FUNCTION("columnMajorAccess-4D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
  for(axom::IndexType l = idxBegin[3]; l < idxEnd[3]; ++l)
  {
    for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
    {
      for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
      {
        for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
        {
          array(i, j, k, l) += 100;
        }
      }
    }
  }
#else
  for(axom::IndexType l = idxBegin[3]; l < idxEnd[3]; ++l)
  {
    for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
    {
      for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
      {
        for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
        {
          array(i, j, k, l) += 100;
        }
      }
    }
  }
#endif
}

//
// Dynamic ordering access tests
//

void runTest_dynamicAccess(axom::ArrayView<double, 1>& array)
{
  AXOM_PERF_MARK_FUNCTION("dynamicAccess-1D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    array[i] += 1000;
  }
#else
  for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
  {
    array[i] += 1000;
  }
#endif
}

void runTest_dynamicAccess(axom::ArrayView<double, 2>& array)
{
  AXOM_PERF_MARK_FUNCTION("dynamicAccess-2D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
  const auto& indexer = array.indexer();
  const auto& slowestDirs = indexer.slowestDirs();
  const axom::StackArray<axom::IndexType, 2> begins {idxBegin[slowestDirs[0]],
                                                     idxBegin[slowestDirs[1]]};
  const axom::StackArray<axom::IndexType, 2> ends {idxEnd[slowestDirs[0]],
                                                   idxEnd[slowestDirs[1]]};
#ifdef AXOM_USE_RAJA
  axom::StackArray<axom::IndexType, 2> idx;
  axom::IndexType& m = idx[slowestDirs[0]];
  axom::IndexType& n = idx[slowestDirs[1]];
  for(m = begins[0]; m < ends[0]; ++m)
  {
    for(n = begins[1]; n < ends[1]; ++n)
    {
      array[idx] += 1000;
    }
  }
#else
  axom::StackArray<axom::IndexType, 2> idx;
  axom::IndexType& m = idx[slowestDirs[0]];
  axom::IndexType& n = idx[slowestDirs[1]];
  for(m = begins[0]; m < ends[0]; ++m)
  {
    for(n = begins[1]; n < ends[1]; ++n)
    {
      array[idx] += 1000;
    }
  }
#endif
}

void runTest_dynamicAccess(axom::ArrayView<double, 3>& array)
{
  AXOM_PERF_MARK_FUNCTION("dynamicAccess-3D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
  const auto& indexer = array.indexer();
  const auto& slowestDirs = indexer.slowestDirs();
  const axom::StackArray<axom::IndexType, 3> begins {idxBegin[slowestDirs[0]],
                                                     idxBegin[slowestDirs[1]],
                                                     idxBegin[slowestDirs[2]]};
  const axom::StackArray<axom::IndexType, 3> ends {idxEnd[slowestDirs[0]],
                                                   idxEnd[slowestDirs[1]],
                                                   idxEnd[slowestDirs[2]]};
#ifdef AXOM_USE_RAJA
  axom::StackArray<axom::IndexType, 3> idx;
  axom::IndexType& m = idx[slowestDirs[0]];
  axom::IndexType& n = idx[slowestDirs[1]];
  axom::IndexType& o = idx[slowestDirs[2]];
  for(m = begins[0]; m < ends[0]; ++m)
  {
    for(n = begins[1]; n < ends[1]; ++n)
    {
      for(o = begins[2]; o < ends[2]; ++o)
      {
        array[idx] += 1000;
      }
    }
  }
#else
  axom::StackArray<axom::IndexType, 3> idx;
  axom::IndexType& m = idx[slowestDirs[0]];
  axom::IndexType& n = idx[slowestDirs[1]];
  axom::IndexType& o = idx[slowestDirs[2]];
  for(m = begins[0]; m < ends[0]; ++m)
  {
    for(n = begins[1]; n < ends[1]; ++n)
    {
      for(o = begins[2]; o < ends[2]; ++o)
      {
        array[idx] += 1000;
      }
    }
  }
#endif
}

void runTest_dynamicAccess(axom::ArrayView<double, 4>& array)
{
  AXOM_PERF_MARK_FUNCTION("dynamicAccess-4D");
  const auto idxBegin = params.idxBegin;
  const auto idxEnd = params.idxEnd;
  const auto& indexer = array.indexer();
  const auto& slowestDirs = indexer.slowestDirs();
  const axom::StackArray<axom::IndexType, 4> begins {idxBegin[slowestDirs[0]],
                                                     idxBegin[slowestDirs[1]],
                                                     idxBegin[slowestDirs[2]],
                                                     idxBegin[slowestDirs[3]]};
  const axom::StackArray<axom::IndexType, 4> ends {idxEnd[slowestDirs[0]],
                                                   idxEnd[slowestDirs[1]],
                                                   idxEnd[slowestDirs[2]],
                                                   idxEnd[slowestDirs[3]]};
#ifdef AXOM_USE_RAJA
  axom::StackArray<axom::IndexType, 4> idx;
  axom::IndexType& m = idx[slowestDirs[0]];
  axom::IndexType& n = idx[slowestDirs[1]];
  axom::IndexType& o = idx[slowestDirs[2]];
  axom::IndexType& p = idx[slowestDirs[3]];
  for(m = begins[0]; m < ends[0]; ++m)
  {
    for(n = begins[1]; n < ends[1]; ++n)
    {
      for(o = begins[2]; o < ends[2]; ++o)
      {
        for(p = begins[3]; p < ends[3]; ++p)
        {
          array[idx] += 1000;
        }
      }
    }
  }
#else
  axom::StackArray<axom::IndexType, 4> idx;
  axom::IndexType& m = idx[slowestDirs[0]];
  axom::IndexType& n = idx[slowestDirs[1]];
  axom::IndexType& o = idx[slowestDirs[2]];
  axom::IndexType& p = idx[slowestDirs[3]];
  for(m = begins[0]; m < ends[0]; ++m)
  {
    for(n = begins[1]; n < ends[1]; ++n)
    {
      for(o = begins[2]; o < ends[2]; ++o)
      {
        for(p = begins[3]; p < ends[3]; ++p)
        {
          array[idx] += 1000;
        }
      }
    }
  }
#endif
}

//
// Make test Array objects.
//

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
void makeArray(axom::Array<double>& ar, axom::ArrayView<double, DIM>& view)
{
  assert(DIM <= params.shape.size());

  axom::StackArray<axom::IndexType, DIM> paddedShape;
  for(int d = 0; d < DIM; ++d)
  {
    paddedShape[d] = params.paddedShape[d];
  }

  ar.resize(params.paddedSize);

  axom::ArrayIndexer<axom::IndexType, DIM> indexer;
  if(!params.slowestDirections.empty())
  {
    axom::StackArray<std::uint16_t, DIM> slowestDirections;
    for(int d = 0; d < DIM; ++d)
    {
      slowestDirections[d] = params.slowestDirections[d];
    }
    indexer.initializeShape(paddedShape, slowestDirections);
  }
  else
  {
    indexer.initializeShape(paddedShape, params.strideOrder);
  }
  view = axom::ArrayView<double, DIM>(ar.data(), paddedShape, indexer.strides());

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

  std::cout << "Real-to-padded size: " << params.realSize << '/'
            << params.paddedSize << " = "
            << double(params.realSize) / params.paddedSize << std::endl;

  auto count = array.size();
  auto* ptr = array.data();
  for(axom::IndexType i = 0; i < count; ++i)
  {
    ptr[i] = double(i * 1000000);
  }

  axom::utilities::Timer sequentialTimer(false);
  sequentialTimer.start();
  runTest_sequentialAccess(array);
  sequentialTimer.stop();
  std::cout << "Sequential time   " << sequentialTimer.elapsedTimeInSec()
            << " seconds, base" << std::endl;

  axom::utilities::Timer rowMajorTimer(false);
  rowMajorTimer.start();
  runTest_rowMajorAccess(array);
  rowMajorTimer.stop();
  std::cout << "Row-major time    " << rowMajorTimer.elapsedTimeInSec()
            << " seconds, " << std::setprecision(3)
            << rowMajorTimer.elapsedTimeInSec() /
      sequentialTimer.elapsedTimeInSec()
            << 'x' << std::endl;

  axom::utilities::Timer columnMajorTimer(false);
  columnMajorTimer.start();
  runTest_columnMajorAccess(array);
  columnMajorTimer.stop();
  std::cout << "Column-major time " << columnMajorTimer.elapsedTimeInSec()
            << " seconds, " << std::setprecision(3)
            << columnMajorTimer.elapsedTimeInSec() /
      sequentialTimer.elapsedTimeInSec()
            << 'x' << std::endl;

  axom::utilities::Timer dynamicTimer(false);
  dynamicTimer.start();
  runTest_dynamicAccess(array);
  dynamicTimer.stop();
  std::cout << "Dynamic time      " << dynamicTimer.elapsedTimeInSec()
            << " seconds, " << std::setprecision(3)
            << dynamicTimer.elapsedTimeInSec() /
      sequentialTimer.elapsedTimeInSec()
            << 'x' << std::endl;
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
