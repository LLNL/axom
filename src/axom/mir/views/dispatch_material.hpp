// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_MATERIAL_HPP_
#define AXOM_MIR_DISPATCH_MATERIAL_HPP_

#include "axom/mir/views/MaterialView.hpp"
#include "axom/mir/views/NodeArrayView.hpp"

#include <conduit/conduit_blueprint.hpp>

namespace axom
{
namespace mir
{
namespace views
{
/**
 * \brief Dispatch a Conduit node containing a matset to a function as the appropriate type of matset view.
 *
 * \tparam FuncType The function/lambda type that will take the matset.
 *
 * \param matset The node that contains the matset.
 * \param func   The function/lambda that will operate on the matset view.
 */
template <typename FuncType>
void dispatch_material(const conduit::Node &matset, FuncType &&func)
{
  constexpr static size_t MaxMaterials = 20;

  if(conduit::blueprint::mesh::matset::is_uni_buffer(matset))
  {
    IndexNode_to_ArrayView_same(
      matset["material_ids"],
      matset["sizes"],
      matset["offsets"],
      matset["indices"],
      [&](auto material_ids, auto sizes, auto offsets, auto indices) {
        FloatNode_to_ArrayView(matset["volume_fractions"], [&](auto volume_fractions) {
          using IndexType = typename decltype(material_ids)::value_type;
          using FloatType = typename decltype(volume_fractions)::value_type;

          UnibufferMaterialView<IndexType, FloatType, MaxMaterials> matsetView;
          matsetView.set(material_ids, volume_fractions, sizes, offsets, indices);
          func(matsetView);
        });
      });
  }
#if 0
  else if(conduit::blueprint::mesh::matset::is_multi_buffer(matset))
  {
    const conduit::Node &volume_fractions = matset.fetch_existing("volume_fractions");
    const conduit::Node &n_firstValues = volume_fractions[0].fetch_existing("values");
    const conduit::Node &n_firstIndices = volume_fractions[0].fetch_existing("indices");
    IndexNode_To_ArrayView(n_firstIndices, [&](auto firstIndices)
    {
      FloatNode_To_ArrayView(n_firstValues, [&](auto firstValues)
      {
        using IntView = decltype(firstIndices);
        using IntElement = typename IntView::value_type;
        using FloatView = decltype(firstValues);
        using FloatElement = typename FloatView::value_type;

        MultiBufferMaterialView<IntElement, FloatElement, MaxMaterials> matsetView;

        for(conduit::index_t i = 0; i < volume_fractions.number_of_children(); i++)
        {
          const conduit::Node &values = volume_fractions[i].fetch_existing("values");
          const conduit::Node &indices = volume_fractions[i].fetch_existing("indices");

          const IntElement *indices_ptr = indices.value();
          const FloatElement *values_ptr = values.value();

          IntView   indices_view(indices_ptr, indices.dtype().number_of_elements());
          FloatView values_view(values_ptr, values.dtype().number_of_elements());
          matsetView.add(indices_view, values_view);
        }

        func(matsetView);
      });
    });
  }
  else if(conduit::blueprint::mesh::matset::is_element_dominant(matset))
  {
    const conduit::Node &volume_fractions = matset.fetch_existing("volume_fractions");
    const conduit::Node &n_firstValues = volume_fractions[0];
    FloatNode_To_ArrayView(n_firstValues, [&](auto firstValues)
    {
      using FloatView = decltype(firstValues);
      using FloatElement = typename FloatView::value_type;

      ElementDominantMaterialView<axom::IndexType, FloatElement, MaxMaterials> matsetView;

      for(conduit::index_t i = 0; i < volume_fractions.number_of_children(); i++)
      {
        const conduit::Node &values = volume_fractions[i];
        const FloatElement *values_ptr = values.value();
        FloatView values_view(values_ptr, values.dtype().number_of_elements());
        matsetView.add(values_view);
      }

      func(matsetView);
    });
  }  
  else if(conduit::blueprint::mesh::matset::is_material_dominant(matset))
  {
    const conduit::Node &volume_fractions = matset.fetch_existing("volume_fractions");
    const conduit::Node &element_ids = matset.fetch_existing("element_ids");   
    const conduit::Node &n_firstValues = volume_fractions[0];
    const conduit::Node &n_firstIndices = element_ids[0];

    IndexNode_To_ArrayView(n_firstIndices, [&](auto firstIndices)
    {
      FloatNode_To_ArrayView(n_firstValues, [&](auto firstValues)
      {
        using IntView = decltype(firstIndices);
        using IntElement = typename IntView::value_type;
        using FloatView = decltype(firstValues);
        using FloatElement = typename FloatView::value_type;

        MaterialDominantMaterialView<IntElement, FloatElement, MaxMaterials> matsetView;

        for(conduit::index_t i = 0; i < volume_fractions.number_of_children(); i++)
        {
          const conduit::Node &indices = element_ids[i];
          const conduit::Node &values = volume_fractions[i];

          const IntElement *indices_ptr = indices.value();
          const FloatElement *values_ptr = values.value();

          IntView   indices_view(indices_ptr, indices.dtype().number_of_elements());
          FloatView values_view(values_ptr, values.dtype().number_of_elements());
          matsetView.add(indices_view, values_view);
        }

        func(matsetView);
      });
    });
  }
#endif
}

}  // end namespace views
}  // end namespace mir
}  // end namespace axom

#endif
