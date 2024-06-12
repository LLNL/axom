// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_CLIPPING_CLIP_TABLE_MANAGER_HPP_
#define AXOM_MIR_CLIPPING_CLIP_TABLE_MANAGER_HPP_

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/mir/clipping/ClipCases.hpp"

namespace axom
{
namespace mir
{
namespace clipping
{

/**
 * \brief This struct contains data for a clipping table.
 *
 * \tparam ContainerType The container for the clipping data.
 * \tparam SPACE The memory space for the clipping data.
 */
template <typename ContainerType, MemorySpace SPACE = MemorySpace::Dynamic>
struct ClipTableBase
{
  using IntContainerType = ContainerType<int, 1, SPACE>;
  using UInt8ContainerType = ContainerType<uint8, 1, SPACE>;

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
  UInt8ContainerType getShape(size_t caseId, size_t shapeId) const
  {
    assert(caseId < m_shapes.size());
    assert(shapeId < shapesForCase(caseId));

    const uint8 *shapeStart = m_table.data() + m_offsets[caseId];
    size_t shapeLen = 0;
    for(int i = 0; i < shapeId; i++)
    {
      shapeLen = advance(*shapeStart);
      shapeStart += shapeLen;
    }
    shapeLen = advance(*shapeStart);
    return UInt8ContainerType(shapeStart, shapeLen);
  }

  /**
   * \brief Given the input shape, return how many values to advance to get to the next shape.
   *
   * \param shape The shape type.
   *
   * \return The number of values to advance.
   */
  AXOM_HOST_DEVICE
  size_t advance(uint8 shape) const
  {
    if(shape == ST_TRI)
      retval = 2 + 3;
    else if(shape == ST_QUAD)
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
  
  IntContainerType   m_shapes;
  IntContainerType   m_offsets;
  UInt8ContainerType m_table;
};

/**
 * \brief This class contains data arrays for the clipping table and can produce a view for the data.
 */
template <MemorySpace SPACE = MemorySpace::Dynamic>
struct ClipTable : public ClipTableBase<axom::Array, SPACE>
{
  using ClipTableView = ClipTableBase<axom::ArrayView, SPACE>

  /**
   * \brief Load clipping data into the arrays, moving data as needed.
   *
   * \param n The number of cases in the clip table.
   * \param shapes The number of shapes produced by clipping cases.
   * \param offsets The offset into the table for each clipping case.
   * \param table The clipping table data.
   * \param tableLen The size of the clipping table data.
   */
  void load(size_t n, const int *shapes, const int *offsets, const uint8 *table, size_t tableLen)
  {
    m_shapes = IntContainerType(shapes, n);
    m_offsets = IntContainerType(offsets, n);
    m_table = UInt8ContainerType(table, tableLen);
  }

  /**
   * \brief Create a view to access the table data.
   *
   * \return A view of the table data.
   */
  ClipTableView view() const
  {
    ClipTableView v;
    v.m_shapes = m_shapes.view();
    v.m_offsets = m_offsets.view();
    v.m_table = m_table.view(); 
    return v;
  }

};

/**
 * \brief Manage several clipping tables.
 */
template <MemorySpace SPACE = MemorySpace::Dynamic>
class ClipTableManager
{
public:
  /**
   * \brief Constructor
   */
  ClipTableManager()
  {
    for(size_t shape = ST_MIN; shape < ST_MAX; shape++)
      m_clipTables[shapeToIndex(shape)] = ClipTable<SPACE>();
  }

  /**
   * \brief Return a reference to the clipping table, which is loaded on demand.
   *
   * \param shape The shape type to be retrieved.
   *
   * \return A reference to the clipping table. 
   */
  const ClipTable<SPACE> &operator[](size_t shape)
  {
    assert(shape < ST_MAX);
    const auto index = shapeToIndex(shape);
    if(m_clipTables[index].size() == 0)
    {
      load(shape);
    }
    return m_clipTables[index];
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
  void load(size_t shape)
  {
    const auto index = shapeToIndex(shape);
    if(shape == ST_TRI)
    {
      m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesTri,
                               axom::mir::clipping::visit::numClipShapesTri,
                               axom::mir::clipping::visit::startClipShapesTri,
                               axom::mir::clipping::visit::clipShapesTri,
                               sizeof(axom::mir::clipping::visit::clipShapesTri) / sizeof(unsigned char));
    }
    else if(shape == ST_QUAD)
    {
      m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesQua,
                               axom::mir::clipping::visit::numClipShapesQua,
                               axom::mir::clipping::visit::startClipShapesQua,
                               axom::mir::clipping::visit::clipShapesQua,
                               sizeof(axom::mir::clipping::visit::clipShapesQua) / sizeof(unsigned char));
    }
    else if(shape == ST_TET)
    {
      m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesTet,
                               axom::mir::clipping::visit::numClipShapesTet,
                               axom::mir::clipping::visit::startClipShapesTet,
                               axom::mir::clipping::visit::clipShapesTet,
                               sizeof(axom::mir::clipping::visit::clipShapesTet) / sizeof(unsigned char));
    }
    else if(shape == ST_PYR)
    {
      m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesPyr,
                               axom::mir::clipping::visit::numClipShapesPyr,
                               axom::mir::clipping::visit::startClipShapesPyr,
                               axom::mir::clipping::visit::clipShapesPyr,
                               sizeof(axom::mir::clipping::visit::clipShapesTet) / sizeof(unsigned char));
    }
    else if(shape == ST_WDG)
    {
      m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesWdg,
                               axom::mir::clipping::visit::numClipShapesWdg,
                               axom::mir::clipping::visit::startClipShapesWdg,
                               axom::mir::clipping::visit::clipShapesWdg,
                               sizeof(axom::mir::clipping::visit::clipShapesWdg) / sizeof(unsigned char));
    }
    else if(shape == ST_HEX)
    {
      m_clipTables[index].load(axom::mir::clipping::visit::numClipCasesHex,
                               axom::mir::clipping::visit::numClipShapesHex,
                               axom::mir::clipping::visit::startClipShapesHex,
                               axom::mir::clipping::visit::clipShapesHex,
                               sizeof(axom::mir::clipping::visit::clipShapesHex) / sizeof(unsigned char));
    }
  }

  ClipTable<SPACE> m_clipTables[ST_MAX - ST_MIN];
};

} // end namespace clipping
} // end namespace mir
} // end namespace axom

#endif
