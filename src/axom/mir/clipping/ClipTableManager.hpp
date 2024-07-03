// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_CLIPPING_CLIP_TABLE_MANAGER_HPP_
#define AXOM_MIR_CLIPPING_CLIP_TABLE_MANAGER_HPP_

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/mir/clipping/ClipCases.h"

namespace axom
{
namespace mir
{
namespace clipping
{

/**
 * \brief This struct contains data for a clipping table.
 *
 * \tparam ClipIndexContainerType The container for int clipping data.
 * \tparam ClipDataContainerType The container for uint8 clipping data.
 */
template <typename ClipIndexContainerType, typename ClipDataContainerType>
struct ClipTableBase
{
  using ClipDataType = typename ClipDataContainerType::value_type;

  // Q: Do I need to explicitly provide a constructor to get it marked as AXOM_HOST_DEVICE?

  /**
   * \brief Return the number of cases for the clipping table.
   *
   * \return The number of cases for the clipping table.
   */
  AXOM_HOST_DEVICE
  size_t size() const { return m_shapes.size(); }

  /**
   * \brief Return the number of shapes for a given clipping case.
   *
   * \param caseId The index of the clipping case.
   *
   * \return The number of shapes in the clipping case.
   */
  AXOM_HOST_DEVICE
  size_t shapesForCase(size_t caseId) const
  {
    assert(caseId < m_shapes.size());
    return m_shapes[caseId];
  }

  /**
   * \brief Return data for the requested shape in the clipping case.
   *
   * \param index The index of the clipping case.
   *
   * \return A container that holds the shape data for the case, probably a view.
   */
  AXOM_HOST_DEVICE
  ClipDataContainerType getShape(size_t caseId, size_t shapeId) const
  {
    assert(caseId < m_shapes.size());
    assert(shapeId < shapesForCase(caseId));

    const auto *shapeStart = m_table.data() + m_offsets[caseId];
    size_t shapeLen = 0;
    for(size_t i = 0; i < shapeId; i++)
    {
      shapeLen = advance(*shapeStart);
      shapeStart += shapeLen;
    }
    shapeLen = advance(*shapeStart);
    return ClipDataContainerType(shapeStart, shapeLen);
  }

  /**
   * \brief Given the input shape, return how many values to advance to get to the next shape.
   *
   * \param shape The shape type.
   *
   * \return The number of values to advance.
   */
  AXOM_HOST_DEVICE
  size_t advance(ClipDataType shape) const
  {
    size_t retval = 0;
    if(shape == ST_TRI)
      retval = 2 + 3;
    else if(shape == ST_QUA)
      retval = 2 + 4;
    else if(shape == ST_TET)
      retval = 2 + 4;
    else if(shape == ST_PYR)
      retval = 2 + 5;
    else if(shape == ST_WDG)
      retval = 2 + 6;
    else if(shape == ST_HEX)
      retval = 2 + 8;
    return retval;
  }
  
  ClipIndexContainerType m_shapes;
  ClipIndexContainerType m_offsets;
  ClipDataContainerType  m_table;
};

/**
 * \brief This class contains data arrays for the clipping table and can produce a view for the data.
 */
template <typename ExecSpace>
struct ClipTable : public ClipTableBase<axom::Array<int>, axom::Array<unsigned char>>
{
  using IntContainerType = axom::Array<int>;
  using Uint8ContainerType = axom::Array<unsigned char>;
  using IntViewType = axom::ArrayView<int>;
  using Uint8ViewType = axom::ArrayView<unsigned char>;

  using SuperClass = ClipTableBase<IntContainerType, Uint8ContainerType>;
  using ClipTableView = ClipTableBase<IntViewType, Uint8ViewType>;

  /**
   * \brief Load clipping data into the arrays, moving data as needed.
   *
   * \param n The number of cases in the clip table.
   * \param shapes The number of shapes produced by clipping cases.
   * \param offsets The offset into the table for each clipping case.
   * \param table The clipping table data.
   * \param tableLen The size of the clipping table data.
   */
  void load(size_t n, const int *shapes, const int *offsets, const unsigned char *table, size_t tableLen)
  {
    const int allocatorID = execution_space<ExecSpace>::allocatorID();

    // Allocate space.
    SuperClass::m_shapes = IntContainerType(n, n, allocatorID);
    SuperClass::m_offsets = IntContainerType(n, n, allocatorID);
    SuperClass::m_table = Uint8ContainerType(tableLen, tableLen, allocatorID);

    // Copy data to the arrays.
    axom::copy(SuperClass::m_shapes.data(), shapes, n * sizeof(int));
    axom::copy(SuperClass::m_offsets.data(), offsets, n * sizeof(int));
    axom::copy(SuperClass::m_table.data(), table, tableLen * sizeof(unsigned char));
  }

  /**
   * \brief Create a view to access the table data.
   *
   * \return A view of the table data.
   */
  ClipTableView view() const
  {
    ClipTableView v;
    v.m_shapes = SuperClass::m_shapes.view();
    v.m_offsets = SuperClass::m_offsets.view();
    v.m_table = SuperClass::m_table.view(); 
    return v;
  }

};

/**
 * \brief Manage several clipping tables.
 */
template <typename ExecSpace>
class ClipTableManager
{
public:
  /**
   * \brief Constructor
   */
  ClipTableManager()
  {
    for(size_t shape = ST_MIN; shape < ST_MAX; shape++)
      m_clipTables[shapeToIndex(shape)] = ClipTable<ExecSpace>();
  }

  /**
   * \brief Return a reference to the clipping table, which is loaded on demand.
   *
   * \param shape The shape type to be retrieved.
   *
   * \return A reference to the clipping table. 
   */
  const ClipTable<ExecSpace> &operator[](size_t shape)
  {
    const auto index = shapeToIndex(shape);
    assert(shape < ST_MAX);
    assert(index >= 0);
    load(shape, 0);
    return m_clipTables[index];
  }

  /**
   * \brief Load tables based on dimension.
   */
  void load(int dim)
  {
    for(const auto shape : shapes(dim))
      load(shape, 0);
  }

  /**
   * \brief Return a vector of clipping shape ids for the given dimension.
   *
   * \param The spatial dimension.
   *
   * \return A vector of clipping shape ids.
   */
  std::vector<size_t> shapes(int dim) const
  {
    std::vector<size_t> s;
    if(dim == 2)
      s = std::vector<size_t>{ST_TRI, ST_QUA};
    else if(dim == 3)
      s = std::vector<size_t>{ST_TET, ST_PYR, ST_WDG, ST_HEX};
    return s;
  }
private:
  /**
   * \brief Turn a shape into an table index.
   *
   * \param shape The shape type ST_XXX.
   *
   * \return An index into the m_clipTables array.
   */
  size_t shapeToIndex(size_t shape) const
  {
    return shape - ST_MIN;
  }

  /**
   * \brief Load the clipping table for a shape.
   *
   * \param shape The shape whose table will be loaded.
   */
  void load(size_t shape, int)
  {
    const auto index = shapeToIndex(shape);
    if(m_clipTables[index].size() == 0)
    {
      if(shape == ST_TRI)
      {
        m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesTri,
                                 axom::mir::clipping::visit::numClipShapesTri,
                                 axom::mir::clipping::visit::startClipShapesTri,
                                 axom::mir::clipping::visit::clipShapesTri,
                                 axom::mir::clipping::visit::clipShapesTriSize);
      }
      else if(shape == ST_QUA)
      {
        m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesQua,
                                 axom::mir::clipping::visit::numClipShapesQua,
                                 axom::mir::clipping::visit::startClipShapesQua,
                                 axom::mir::clipping::visit::clipShapesQua,
                                 axom::mir::clipping::visit::clipShapesQuaSize);
      }
      else if(shape == ST_TET)
      {
        m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesTet,
                                 axom::mir::clipping::visit::numClipShapesTet,
                                 axom::mir::clipping::visit::startClipShapesTet,
                                 axom::mir::clipping::visit::clipShapesTet,
                                 axom::mir::clipping::visit::clipShapesTetSize);
      }
      else if(shape == ST_PYR)
      {
        m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesPyr,
                                 axom::mir::clipping::visit::numClipShapesPyr,
                                 axom::mir::clipping::visit::startClipShapesPyr,
                                 axom::mir::clipping::visit::clipShapesPyr,
                                 axom::mir::clipping::visit::clipShapesTetSize);
      }
      else if(shape == ST_WDG)
      {
        m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesWdg,
                                 axom::mir::clipping::visit::numClipShapesWdg,
                                 axom::mir::clipping::visit::startClipShapesWdg,
                                 axom::mir::clipping::visit::clipShapesWdg,
                                 axom::mir::clipping::visit::clipShapesWdgSize);
      }
      else if(shape == ST_HEX)
      {
        m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesHex,
                                 axom::mir::clipping::visit::numClipShapesHex,
                                 axom::mir::clipping::visit::startClipShapesHex,
                                 axom::mir::clipping::visit::clipShapesHex,
                                 axom::mir::clipping::visit::clipShapesHexSize);
      }
    }
  }

  ClipTable<ExecSpace> m_clipTables[ST_MAX - ST_MIN];
};

} // end namespace clipping
} // end namespace mir
} // end namespace axom

#endif
