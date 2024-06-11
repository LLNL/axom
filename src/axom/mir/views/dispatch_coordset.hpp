// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_DISPATCH_COORDSET_HPP_
#define AXOM_MIR_DISPATCH_COORDSET_HPP_

#include "axom/mir/views/ExplicitCoordsetView.hpp"
#include "axom/mir/views/UniformCoordsetView.hpp"
#include "axom/mir/views/RectilinearCoordsetView.hpp"
#include "axom/mir/views/NodeArrayView.hpp"

namespace axom
{
namespace mir
{
namespace views
{

/**
 * \brief Given a Conduit/Blueprint coordset, create an appropriate view and
 *        call the supplied function, passing the coordset view to it.
 *
 * \tparam FuncType The type of the function / lambda to invoke. It is expected
 *                  that the callable accepts an auto argument for a coordset view.
 *
 * \param coordset The Conduit node that contains the coordset.
 * \param func     The function/lambda to invoke using the coordset view.
 */
template <typename FuncType>
void dispatch_coordset(const conduit::Node &coordset, FuncType &&func)
{
  const std::string cstype = coordset["type"].as_string();
  if(cstype == "uniform")
  {   
    const conduit::Node &n_dims = coordset["dims"];
    const conduit::index_t ndims = n_dims.dtype().number_of_elements();
    if(ndims == 2)
    {
      axom::StackArray<axom::IndexType, 2> dims;
      axom::StackArray<double, 2> origin{0., 0.}, spacing{1., 1.};
      for(int i = 0; i < ndims; i++)
      {
        dims[i] = n_dims.as_int_accessor()[i];
        if(coordset.has_child("origin"))
          origin[i] = coordset["origin"].as_double_accessor()[i];
        if(coordset.has_child("spacing"))
          spacing[i] = coordset["spacing"].as_double_accessor()[i];
      }

      UniformCoordsetView<double, 2> coordView(dims, origin, spacing);
      func(coordView);
    }
    else if(ndims == 3)
    {
      axom::StackArray<axom::IndexType, 3> dims;
      axom::StackArray<double, 3> origin{0., 0., 0.}, spacing{1., 1., 1.};
      for(int i = 0; i < ndims; i++)
      {
        dims[i] = n_dims.as_int_accessor()[i];
        if(coordset.has_child("origin"))
          origin[i] = coordset["origin"].as_double_accessor()[i];
        if(coordset.has_child("spacing"))
          spacing[i] = coordset["spacing"].as_double_accessor()[i];
      }

      UniformCoordsetView<double, 3> coordView(dims, origin, spacing);
      func(coordView);
    }   
  }
  else if(cstype == "rectilinear")
  {
    const conduit::Node &values = coordset["values"];
    if(values.number_of_children() == 2)
    {
      axom::mir::views::FloatNode_to_ArrayView_same(values[0], values[1], [&](auto xView, auto yView)
      {
        RectilinearCoordsetView2<typename decltype(xView)::value_type> coordView(xView, yView);
        func(coordView);
      });
    }
    else if(values.number_of_children() == 3)
    {
      axom::mir::views::FloatNode_to_ArrayView_same(values[0], values[1], values[2], [&](auto xView, auto yView, auto zView)
      {
        RectilinearCoordsetView3<typename decltype(xView)::value_type> coordView(xView, yView, zView);
        func(coordView);
      });
    }
  }
  else if(cstype == "explicit")
  {
    // TODO: get the axis names.
    // TODO: support strided/structured.
    const conduit::Node &values = coordset["values"];
    if(values.number_of_children() == 2)
    {
      axom::mir::views::FloatNode_to_ArrayView_same(values[0], values[1], [&](auto xView, auto yView)
      {
        ExplicitCoordsetView2<typename decltype(xView)::value_type> coordView(xView, yView);
        func(coordView);
      });
    }
    else if(values.number_of_children() == 3)
    {
      axom::mir::views::FloatNode_to_ArrayView_same(values[0], values[1], values[2], [&](auto xView, auto yView, auto zView)
      {
        ExplicitCoordsetView3<typename decltype(xView)::value_type> coordView(xView, yView, zView);
        func(coordView);
      });
    }
  }
}

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
