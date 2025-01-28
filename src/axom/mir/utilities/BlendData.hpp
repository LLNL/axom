// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_BLEND_DATA_HPP_
#define AXOM_MIR_BLEND_DATA_HPP_

#include "axom/core.hpp"

namespace axom
{
namespace mir
{
namespace utilities
{
namespace blueprint
{
/*!
 * \brief This class contains views of blend data. Blend data lets is make new
 *        nodal fields and coordsets. The field data are sampled using m_originalIdsView
 *        which is a compact list of the original node ids that we want to preserve
 *        without any blending. This stream is followed by a second stream of data
 *        made using the field and the blend groups. Each blend group has
 *        m_blendGroupSizesView[i] elements, starts at m_blendGroupStartView[i] and 
 *        uses values from the m_blendIdsView, m_blendCoeffView members to blend the
 *        data values.
 *
 */
struct BlendData
{
  using IndexView = axom::ArrayView<IndexType>;
  using CoeffView = axom::ArrayView<float>;

  IndexView m_originalIdsView {};      //!< Indices of original node ids.
  IndexView m_selectedIndicesView {};  //!< Indices of the selected blend groups.
  IndexView m_blendGroupSizesView {};  //!< Number of ids/weights in each blend group.
  IndexView m_blendGroupStartView {};  //!< Starting offset for a blend group in ids.
  IndexView m_blendIdsView {};         //!< Ids that make up blend groups.
  CoeffView m_blendCoeffView {};       //!< Weights that make up blend groups.
};

}  // end namespace blueprint
}  // end namespace utilities
}  // end namespace mir
}  // end namespace axom

#endif
