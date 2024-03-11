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

  RuntimePolicy runtimePolicy = RuntimePolicy::seq;

private:
  bool _verboseOutput {false};
  const std::map<std::string, int> strideValidator {
    {"row", int(axom::ArrayStrideOrder::ROW)},
    {"col", int(axom::ArrayStrideOrder::COLUMN)}};

public:
  bool isVerbose() const { return _verboseOutput; }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-p, --policy", runtimePolicy)
      ->description("Set runtime policy for test")
      ->capture_default_str()
      ->transform(
        axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

    app.add_flag("-v,--verbose,!--no-verbose", _verboseOutput)
      ->description("Enable/disable verbose output")
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

template <int DIM, typename ExecSpace>
class ArrayIndexerPerfTester
{
public:
  using Element_t = std::uint64_t;
  const Element_t m_baseFactor = 1000000;
  Element_t m_testAccumulation = 0;
  const Element_t sequentiaTestAdd = 1;
  const Element_t rowTestAdd = 10;
  const Element_t columnTestAdd = 100;
  const Element_t dynamicTestAdd = 1000;

  /*!
  @brief Time the sequential access of every element of an array.

  This is the fastest we expect to visit every element.
*/
  void runTest_sequentialAccess(axom::ArrayView<Element_t, DIM>& array)
  {
    auto testAdd = sequentiaTestAdd;
    m_testAccumulation += testAdd;
    auto count = array.size();
    auto* ptr = array.data();
#ifdef AXOM_USE_RAJA
    axom::for_all<ExecSpace>(
      0,
      count,
      AXOM_LAMBDA(axom::IndexType i) { ptr[i] += testAdd; });
#else
    for(axom::IndexType i = 0; i < count; ++i)
    {
      ptr[i] += testAdd;
    }
#endif
  }

  /*!
  Methods to time the access of every element of an array.

  Multidimensional loops are capable of skipping ghost
  layers, but the sequential loop used for the baseline
  performance doesn't have logic to skip them.
*/

  //
  // Row-major access tests
  //

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 1>::type runTest_rowMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("rowMajorAccess-1D");
    auto testAdd = rowTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
    axom::for_all<ExecSpace>(
      idxBegin[0],
      idxEnd[0],
      AXOM_LAMBDA(axom::IndexType i) { array[i] += testAdd; });
#else
    for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
    {
      array[i] += testAdd;
    }
#endif
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type runTest_rowMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("rowMajorAccess-2D");
    auto testAdd = rowTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
    using EXEC_POL =
      typename axom::internal::nested_for_exec<ExecSpace>::loop2d_policy;
    RAJA::RangeSegment iRange(idxBegin[0], idxEnd[0]);
    RAJA::RangeSegment jRange(idxBegin[1], idxEnd[1]);
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(jRange, iRange),
      AXOM_LAMBDA(axom::IndexType j, axom::IndexType i) {
        array(i, j) += testAdd;
      });
#else
    for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
    {
      for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
      {
        array(i, j) += testAdd;
      }
    }
#endif
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type runTest_rowMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("rowMajorAccess-3D");
    auto testAdd = rowTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
    using EXEC_POL =
      typename axom::internal::nested_for_exec<ExecSpace>::loop3d_policy;
    RAJA::RangeSegment iRange(idxBegin[0], idxEnd[0]);
    RAJA::RangeSegment jRange(idxBegin[1], idxEnd[1]);
    RAJA::RangeSegment kRange(idxBegin[2], idxEnd[2]);
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(kRange, jRange, iRange),
      AXOM_LAMBDA(axom::IndexType k, axom::IndexType j, axom::IndexType i) {
        array(i, j, k) += testAdd;
      });
#else
    for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
    {
      for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
      {
        for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
        {
          array(i, j, k) += testAdd;
        }
      }
    }
#endif
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 4>::type runTest_rowMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("rowMajorAccess-4D");
    auto testAdd = rowTestAdd;
    m_testAccumulation += testAdd;

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
            array(i, j, k, l) += testAdd;
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
            array(i, j, k, l) += testAdd;
          }
        }
      }
    }
#endif
  }

  //
  // Colunn-major access tests
  //

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 1>::type runTest_columnMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("columnMajorAccess-1D");
    auto testAdd = columnTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
    axom::for_all<ExecSpace>(
      idxBegin[0],
      idxEnd[0],
      AXOM_LAMBDA(axom::IndexType i) { array[i] += testAdd; });
#else
    for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
    {
      array[i] += testAdd;
    }
#endif
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type runTest_columnMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("columnMajorAccess-2D");
    auto testAdd = columnTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
    using EXEC_POL =
      typename axom::internal::nested_for_exec<ExecSpace>::loop2d_policy;
    RAJA::RangeSegment iRange(idxBegin[0], idxEnd[0]);
    RAJA::RangeSegment jRange(idxBegin[1], idxEnd[1]);
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(iRange, jRange),
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j) {
        array(i, j) += testAdd;
      });
#else
    for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
    {
      for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
      {
        array(i, j) += testAdd;
      }
    }
#endif
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type runTest_columnMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("columnMajorAccess-3D");
    auto testAdd = columnTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
    using EXEC_POL =
      typename axom::internal::nested_for_exec<ExecSpace>::loop3d_policy;
    RAJA::RangeSegment iRange(idxBegin[0], idxEnd[0]);
    RAJA::RangeSegment jRange(idxBegin[1], idxEnd[1]);
    RAJA::RangeSegment kRange(idxBegin[2], idxEnd[2]);
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(iRange, jRange, kRange),
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j, axom::IndexType k) {
        array(i, j, k) += testAdd;
      });
#else
    for(axom::IndexType k = idxBegin[2]; k < idxEnd[2]; ++k)
    {
      for(axom::IndexType j = idxBegin[1]; j < idxEnd[1]; ++j)
      {
        for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
        {
          array(i, j, k) += testAdd;
        }
      }
    }
#endif
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 4>::type runTest_columnMajorAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("columnMajorAccess-4D");
    auto testAdd = columnTestAdd;
    m_testAccumulation += testAdd;

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
            array(i, j, k, l) += testAdd;
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
            array(i, j, k, l) += testAdd;
          }
        }
      }
    }
#endif
  }

  //
  // Dynamic ordering access tests

  // The dynamic order should match the optimal order,
  // Any performance difference is due to overhead of dynamic
  // nesting of the loops.
  //

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 1>::type runTest_dynamicAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("dynamicAccess-1D");
    auto testAdd = dynamicTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
#ifdef AXOM_USE_RAJA
    axom::for_all<ExecSpace>(
      idxBegin[0],
      idxEnd[0],
      AXOM_LAMBDA(axom::IndexType i) { array[i] += testAdd; });
#else
    for(axom::IndexType i = idxBegin[0]; i < idxEnd[0]; ++i)
    {
      array[i] += testAdd;
    }
#endif
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type runTest_dynamicAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("dynamicAccess-2D");
    auto testAdd = dynamicTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    const auto& indexer = array.indexer();
    const auto& slowestDirs = indexer.slowestDirs();
    axom::StackArray<axom::IndexType, DIM> invSlowestDirs;
    for(int i = 0; i < DIM; ++i) invSlowestDirs[slowestDirs[i]] = i;
    const axom::StackArray<axom::IndexType, DIM> begins {
      idxBegin[slowestDirs[0]],
      idxBegin[slowestDirs[1]]};
    const axom::StackArray<axom::IndexType, DIM> ends {idxEnd[slowestDirs[0]],
                                                       idxEnd[slowestDirs[1]]};
#ifdef AXOM_USE_RAJA
    using EXEC_POL =
      typename axom::internal::nested_for_exec<ExecSpace>::loop2d_policy;
    RAJA::RangeSegment mRange(begins[0], ends[0]);
    RAJA::RangeSegment nRange(begins[1], ends[1]);
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(nRange, mRange),
      AXOM_LAMBDA(axom::IndexType n, axom::IndexType m) {
        axom::StackArray<axom::IndexType, DIM> idx {m, n};
        auto i = idx[invSlowestDirs[0]];
        auto j = idx[invSlowestDirs[1]];
        array(i, j) += testAdd;
      });
#else
    axom::StackArray<axom::IndexType, DIM> idx;
    axom::IndexType& m = idx[slowestDirs[0]];
    axom::IndexType& n = idx[slowestDirs[1]];
    for(m = begins[0]; m < ends[0]; ++m)
    {
      for(n = begins[1]; n < ends[1]; ++n)
      {
        array[idx] += testAdd;
      }
    }
#endif
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type runTest_dynamicAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("dynamicAccess-3D");
    auto testAdd = dynamicTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    const auto& indexer = array.indexer();
    const auto& slowestDirs = indexer.slowestDirs();
    axom::StackArray<axom::IndexType, DIM> invSlowestDirs;
    for(int i = 0; i < DIM; ++i) invSlowestDirs[slowestDirs[i]] = i;
    const axom::StackArray<axom::IndexType, DIM> begins {
      idxBegin[slowestDirs[0]],
      idxBegin[slowestDirs[1]],
      idxBegin[slowestDirs[2]]};
    const axom::StackArray<axom::IndexType, DIM> ends {idxEnd[slowestDirs[0]],
                                                       idxEnd[slowestDirs[1]],
                                                       idxEnd[slowestDirs[2]]};
#ifdef AXOM_USE_RAJA
    using EXEC_POL =
      typename axom::internal::nested_for_exec<ExecSpace>::loop3d_policy;
    RAJA::RangeSegment mRange(begins[0], ends[0]);
    RAJA::RangeSegment nRange(begins[1], ends[1]);
    RAJA::RangeSegment oRange(begins[2], ends[2]);
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(oRange, nRange, mRange),
      AXOM_LAMBDA(axom::IndexType o, axom::IndexType n, axom::IndexType m) {
        axom::StackArray<axom::IndexType, DIM> idx {m, n, o};
        auto i = idx[invSlowestDirs[0]];
        auto j = idx[invSlowestDirs[1]];
        auto k = idx[invSlowestDirs[2]];
        array(i, j, k) += testAdd;
      });
#else
    axom::StackArray<axom::IndexType, DIM> idx;
    axom::IndexType& m = idx[slowestDirs[0]];
    axom::IndexType& n = idx[slowestDirs[1]];
    axom::IndexType& o = idx[slowestDirs[2]];
    for(m = begins[0]; m < ends[0]; ++m)
    {
      for(n = begins[1]; n < ends[1]; ++n)
      {
        for(o = begins[2]; o < ends[2]; ++o)
        {
          array[idx] += testAdd;
        }
      }
    }
#endif
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 4>::type runTest_dynamicAccess(
    axom::ArrayView<Element_t, DIM>& array)
  {
    AXOM_PERF_MARK_FUNCTION("dynamicAccess-4D");
    auto testAdd = dynamicTestAdd;
    m_testAccumulation += testAdd;

    const auto idxBegin = params.idxBegin;
    const auto idxEnd = params.idxEnd;
    const auto& indexer = array.indexer();
    const auto& slowestDirs = indexer.slowestDirs();
    const axom::StackArray<axom::IndexType, DIM> begins {
      idxBegin[slowestDirs[0]],
      idxBegin[slowestDirs[1]],
      idxBegin[slowestDirs[2]],
      idxBegin[slowestDirs[3]]};
    const axom::StackArray<axom::IndexType, DIM> ends {idxEnd[slowestDirs[0]],
                                                       idxEnd[slowestDirs[1]],
                                                       idxEnd[slowestDirs[2]],
                                                       idxEnd[slowestDirs[3]]};
#ifdef AXOM_USE_RAJA
    axom::StackArray<axom::IndexType, DIM> idx;
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
            array[idx] += testAdd;
          }
        }
      }
    }
#else
    axom::StackArray<axom::IndexType, DIM> idx;
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
            array[idx] += testAdd;
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
  axom::Array<Element_t, DIM> makeArray()
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

    return axom::Array<Element_t, DIM>(shape);
  }

  /*!
  @brief Return an array for testing, dimension DIM,
  sized and ordered according to params values.

  This method allocates a 1D array and puts a muldimensional
  view on the data.  The view supports arbitrary ordering.
*/
  void makeArray(axom::Array<Element_t>& ar, axom::ArrayView<Element_t, DIM>& view)
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
    view =
      axom::ArrayView<Element_t, DIM>(ar.data(), paddedShape, indexer.strides());

    return;
  }

  /*!
  @brief Run test with array shape of the first DIM values in
  params.shape.
*/
  void runTest_dim()
  {
    // Use ArrayView to test, because Array doesn't support
    // arbitrary ordering (yet).
#if 0
  axom::Array<Element_t, DIM> arrayMd = makeArray();
  axom::ArrayView<Element_t, DIM> array = arrayMd.view();
#else
    axom::Array<Element_t> array1D;
    axom::ArrayView<Element_t, DIM> array;
    makeArray(array1D, array);
#endif

    std::cout << "Real-to-padded size: " << params.realSize << '/'
              << params.paddedSize << " = "
              << double(params.realSize) / params.paddedSize << std::endl;

    auto count = array.size();
    auto* ptr = array.data();
    for(axom::IndexType i = 0; i < count; ++i)
    {
      ptr[i] = Element_t(i * m_baseFactor);
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

    // Verify that the elements are touched the correct number of times.
    axom::IndexType matchCount = 0;
    for(axom::IndexType i = 0; i < count; ++i)
    {
      matchCount += (ptr[i] == Element_t(i * m_baseFactor) + m_testAccumulation);
    }
    if(matchCount != params.realSize)
    {
      std::cerr << "Unexpected error in tests: counted match (" << matchCount
                << ") != expected (" << params.realSize << ")" << std::endl;
    }
  }

};  // class ArrayIndexerPerfTester

/*!
  @brief Run test based on dimension specified in params.
*/
template <typename ExecSpace>
void runTest()
{
  if(params.shape.size() == 1)
  {
    ArrayIndexerPerfTester<1, ExecSpace> tester1D;
    tester1D.runTest_dim();
  }
  else if(params.shape.size() == 2)
  {
    ArrayIndexerPerfTester<2, ExecSpace> tester2D;
    tester2D.runTest_dim();
  }
  else if(params.shape.size() == 3)
  {
    ArrayIndexerPerfTester<3, ExecSpace> tester3D;
    tester3D.runTest_dim();
  }
  else if(params.shape.size() == 4)
  {
    ArrayIndexerPerfTester<4, ExecSpace> tester4D;
    tester4D.runTest_dim();
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

  using RuntimePolicy = axom::runtime_policy::Policy;

  std::cout << "Testing runtimePolicy "
            << axom::runtime_policy::policyToName(params.runtimePolicy)
            << std::endl;

  if(params.runtimePolicy == RuntimePolicy::seq)
  {
    runTest<axom::SEQ_EXEC>();
  }
#ifdef AXOM_RUNTIME_POLICY_USE_OPENMP
  else if(params.runtimePolicy == RuntimePolicy::omp)
  {
    runTest<axom::OMP_EXEC>();
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_CUDA
  else if(params.runtimePolicy == RuntimePolicy::cuda)
  {
    runTest<axom::CUDA_EXEC<256>>();
  }
#endif
#ifdef AXOM_RUNTIME_POLICY_USE_HIP
  else if(params.runtimePolicy == RuntimePolicy::hip)
  {
    runTest<axom::HIP_EXEC<256>>();
  }
#endif
  return 0;
}
