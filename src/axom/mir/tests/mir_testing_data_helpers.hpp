// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_TESTING_DATA_HELPERS_HPP_
#define AXOM_MIR_TESTING_DATA_HELPERS_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include <conduit/conduit.hpp>
#include <string>
#include <vector>

namespace axom
{
namespace mir
{
namespace testing
{
namespace data
{

//------------------------------------------------------------------------------

void add_distance(conduit::Node &mesh, float dist = 6.5f)
{
  // Make a new distance field.
  const conduit::Node &n_coordset = mesh["coordsets"][0];
  axom::mir::views::dispatch_coordset(n_coordset, [&](auto coordsetView) {
    mesh["fields/distance/topology"] = "mesh";
    mesh["fields/distance/association"] = "vertex";
    conduit::Node &n_values = mesh["fields/distance/values"];
    const auto nnodes = coordsetView.size();
    n_values.set(conduit::DataType::float32(nnodes));
    float *valuesPtr = static_cast<float *>(n_values.data_ptr());
    for(int index = 0; index < nnodes; index++)
    {
      const auto pt = coordsetView[index];
      float norm2 = 0.f;
      for(int i = 0; i < pt.DIMENSION; i++) norm2 += pt[i] * pt[i];
      valuesPtr[index] = sqrt(norm2) - dist;
    }
  });
}

template <typename Dimensions>
void braid(const std::string &type, const Dimensions &dims, conduit::Node &mesh)
{
  int d[3] = {0, 0, 0};
  for(int i = 0; i < dims.size(); i++) d[i] = dims[i];
  conduit::blueprint::mesh::examples::braid(type, d[0], d[1], d[2], mesh);

  add_distance(mesh);
}

void mixed3d(conduit::Node &mesh)
{
  // clang-format off
  const std::vector<int> conn{{
    // tets
    0,6,1,3,
    3,6,1,9,
    6,7,1,9,
    3,9,1,4,
    9,7,1,4,
    9,7,4,10,
    // pyramids
    1,7,8,2,4,
    11,8,7,10,4,
    2,8,11,5,4,
    // wedges
    6,7,9,12,13,15,
    9,7,10,15,13,16,
    // hex
    7,13,14,8,10,16,17,11
  }};
  const std::vector<int> shapes{{
    0,0,0,0,0,0,
    1,1,1,
    2,2,
    3
  }};
  const std::vector<int> sizes{{
    4,4,4,4,4,4,
    5,5,5,
    6,6,
    8
  }};
  const std::vector<int> offsets{{
    0,4,8,12,16,20,
    24,29,34,
    39,45,
    51
  }};
  constexpr float LOW = -10.f;
  constexpr float MID = 0.f;
  constexpr float HIGH = 10.f;
  const std::vector<float> x{{
    LOW, MID, HIGH, LOW, MID, HIGH, LOW, MID, HIGH, LOW, MID, HIGH, LOW, MID, HIGH, LOW, MID, HIGH
  }};
  const std::vector<float> y{{
    LOW, LOW, LOW, MID, MID, MID, LOW, LOW, LOW, MID, MID, MID, LOW, LOW, LOW, MID, MID, MID
  }};
  const std::vector<float> z{{
    LOW, LOW, LOW, LOW, LOW, LOW, MID, MID, MID, MID, MID, MID, HIGH, HIGH, HIGH, HIGH, HIGH, HIGH
  }};
  // clang-format off

  mesh["coordsets/coords/type"] = "explicit";
  mesh["coordsets/coords/values/x"].set(x);
  mesh["coordsets/coords/values/y"].set(y);
  mesh["coordsets/coords/values/z"].set(z);
  mesh["topologies/mesh/type"] = "unstructured";
  mesh["topologies/mesh/coordset"] = "coords";
  mesh["topologies/mesh/elements/shape"] = "mixed";
  mesh["topologies/mesh/elements/connectivity"].set(conn);
  mesh["topologies/mesh/elements/shapes"].set(shapes);
  mesh["topologies/mesh/elements/sizes"].set(sizes);
  mesh["topologies/mesh/elements/offsets"].set(offsets);
  mesh["topologies/mesh/elements/shape_map/tet"] = 0;
  mesh["topologies/mesh/elements/shape_map/pyramid"] = 1;
  mesh["topologies/mesh/elements/shape_map/wedge"] = 2;
  mesh["topologies/mesh/elements/shape_map/hex"] = 3;

  add_distance(mesh, 0.f);
}

void make_one_hex(conduit::Node &hostMesh)
{
  hostMesh["coordsets/coords/type"] = "explicit";
  hostMesh["coordsets/coords/values/x"].set(
    std::vector<float> {{0., 1., 1., 0., 0., 1., 1., 0.}});
  hostMesh["coordsets/coords/values/y"].set(
    std::vector<float> {{0., 0., 1., 1., 0., 0., 1., 1.}});
  hostMesh["coordsets/coords/values/z"].set(
    std::vector<float> {{0., 0., 0., 0., 1., 1., 1., 1.}});
  hostMesh["topologies/topo/type"] = "unstructured";
  hostMesh["topologies/topo/coordset"] = "coords";
  hostMesh["topologies/topo/elements/shape"] = "hex";
  hostMesh["topologies/topo/elements/connectivity"].set(
    std::vector<int> {{0, 1, 2, 3, 4, 5, 6, 7}});
  hostMesh["topologies/topo/elements/sizes"].set(std::vector<int> {8});
  hostMesh["topologies/topo/elements/offsets"].set(std::vector<int> {0});
  hostMesh["fields/distance/topology"] = "topo";
  hostMesh["fields/distance/association"] = "vertex";
  hostMesh["fields/distance/values"].set(
    std::vector<float> {{1., -1., -1., -1., -1., -1., -1., -1.}});
}

void make_one_tet(conduit::Node &hostMesh)
{
  hostMesh["coordsets/coords/type"] = "explicit";
  hostMesh["coordsets/coords/values/x"].set(std::vector<float> {{0., 0., 1., 0.}});
  hostMesh["coordsets/coords/values/y"].set(std::vector<float> {{0., 0., 0., 1.}});
  hostMesh["coordsets/coords/values/z"].set(std::vector<float> {{0., 1., 0., 0.}});
  hostMesh["topologies/topo/type"] = "unstructured";
  hostMesh["topologies/topo/coordset"] = "coords";
  hostMesh["topologies/topo/elements/shape"] = "tet";
  hostMesh["topologies/topo/elements/connectivity"].set(
    std::vector<int> {{0, 1, 2, 3}});
  hostMesh["topologies/topo/elements/sizes"].set(std::vector<int> {4});
  hostMesh["topologies/topo/elements/offsets"].set(std::vector<int> {0});
  hostMesh["fields/distance/topology"] = "topo";
  hostMesh["fields/distance/association"] = "vertex";
  hostMesh["fields/distance/values"].set(std::vector<float> {{-1., -1., -1., 1.}});
}

void make_one_pyr(conduit::Node &hostMesh)
{
  hostMesh["coordsets/coords/type"] = "explicit";
  hostMesh["coordsets/coords/values/x"].set(
    std::vector<float> {{0., 0., 1., 1., 0.5}});
  hostMesh["coordsets/coords/values/y"].set(
    std::vector<float> {{0., 0., 0., 0., 1.}});
  hostMesh["coordsets/coords/values/z"].set(
    std::vector<float> {{0., 1., 1., 0., 0.5}});
  hostMesh["topologies/topo/type"] = "unstructured";
  hostMesh["topologies/topo/coordset"] = "coords";
  hostMesh["topologies/topo/elements/shape"] = "pyramid";
  hostMesh["topologies/topo/elements/connectivity"].set(
    std::vector<int> {{0, 1, 2, 3, 4}});
  hostMesh["topologies/topo/elements/sizes"].set(std::vector<int> {5});
  hostMesh["topologies/topo/elements/offsets"].set(std::vector<int> {0});
  hostMesh["fields/distance/topology"] = "topo";
  hostMesh["fields/distance/association"] = "vertex";
  hostMesh["fields/distance/values"].set(
    std::vector<float> {{1., 1., -1., -1., -1.}});
}

void make_one_wdg(conduit::Node &hostMesh)
{
  hostMesh["coordsets/coords/type"] = "explicit";
  hostMesh["coordsets/coords/values/x"].set(
    std::vector<float> {{0., 0., 1., 0., 0., 1}});
  hostMesh["coordsets/coords/values/y"].set(
    std::vector<float> {{0., 0., 0., 1., 1., 1.}});
  hostMesh["coordsets/coords/values/z"].set(
    std::vector<float> {{0., 1., 0., 0., 1., 0.}});
  hostMesh["topologies/topo/type"] = "unstructured";
  hostMesh["topologies/topo/coordset"] = "coords";
  hostMesh["topologies/topo/elements/shape"] = "wedge";
  hostMesh["topologies/topo/elements/connectivity"].set(
    std::vector<int> {{0, 1, 2, 3, 4, 5}});
  hostMesh["topologies/topo/elements/sizes"].set(std::vector<int> {6});
  hostMesh["topologies/topo/elements/offsets"].set(std::vector<int> {0});
  hostMesh["fields/distance/topology"] = "topo";
  hostMesh["fields/distance/association"] = "vertex";
  hostMesh["fields/distance/values"].set(
    std::vector<float> {{1., 1., -1., -1., -1., -1.}});
}

template <int NDIMS = 2>
void strided_structured(conduit::Node &hostMesh)
{
  // Total padding in each dimension 
  const conduit::index_t total_elt_pad = 4; // two on each end
  const conduit::index_t total_pt_pad = 3;  // two on the low end, one on the high end

  // Size of the window in the larger mesh
  const conduit::index_t npts_x = 4;
  const conduit::index_t npts_y = 3;
  const conduit::index_t npts_z = (NDIMS == 3) ? 3 : 0;
  const conduit::index_t spz = (NDIMS == 3) ? (npts_z + total_pt_pad) : 0;

  const conduit::index_t nelts_x = npts_x - 1;
  const conduit::index_t nelts_y = npts_y - 1;
  const conduit::index_t nelts_z = (NDIMS == 3) ? (npts_z - 1) : 0;
  const conduit::index_t sez = (NDIMS == 3) ? (nelts_z + total_elt_pad) : 0;

  // Origin: where the data starts in the arrays 
  const conduit::index_t origin_x = 2;
  const conduit::index_t origin_y = 2;
  const conduit::index_t origin_z = (NDIMS == 3) ? 2 : 0;

  conduit::Node desc;
  desc["vertex_data/shape"].set(std::vector<conduit::index_t>{{npts_x + total_pt_pad, npts_y + total_pt_pad, spz}});
  desc["vertex_data/origin"].set(std::vector<conduit::index_t>{{origin_x, origin_y, origin_z}});
  desc["element_data/shape"].set(std::vector<conduit::index_t>{{nelts_x + total_elt_pad, nelts_y + total_elt_pad, sez}});
  desc["element_data/origin"].set(std::vector<conduit::index_t>{{origin_x, origin_y, origin_z}});
  conduit::blueprint::mesh::examples::strided_structured(desc, npts_x, npts_y, npts_z, hostMesh);
}

} // end namespace data
} // end namespace testing
} // end namespace mir
} // end namespace axom

#endif
