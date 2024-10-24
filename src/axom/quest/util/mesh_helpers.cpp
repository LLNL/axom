// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "mesh_helpers.hpp"
#include <axom/sidre.hpp>
#include <axom/slic.hpp>
#include <conduit/conduit_blueprint_mesh.hpp>
#include <iostream>
#include "axom/mint/mesh/UnstructuredMesh.hpp"

namespace axom
{
namespace quest
{
namespace util
{
#ifdef AXOM_USE_MFEM
mfem::Mesh* make_cartesian_mfem_mesh_2D(const primal::BoundingBox<double, 2>& bbox,
                                        const primal::NumericArray<int, 2>& res,
                                        int polynomial_order,
                                        bool reorder_space_filling)
{
  constexpr int DIM = 2;
  const auto range = bbox.range();

  auto* mesh =
    new mfem::Mesh(mfem::Mesh::MakeCartesian2D(res[0],
                                               res[1],
                                               mfem::Element::QUADRILATERAL,
                                               true,
                                               range[0],
                                               range[1],
                                               reorder_space_filling));

  // Offset the mesh to lie w/in the bounding box
  const int NV = mesh->GetNV();
  for(int i = 0; i < NV; ++i)
  {
    double* v = mesh->GetVertex(i);
    for(int d = 0; d < DIM; ++d)
    {
      v[d] += bbox.getMin()[d];
    }
  }

  // Ensure that mesh has high order nodes
  mesh->SetCurvature(polynomial_order);

  return mesh;
}

mfem::Mesh* make_cartesian_mfem_mesh_3D(const primal::BoundingBox<double, 3>& bbox,
                                        const primal::NumericArray<int, 3>& res,
                                        int polynomial_order,
                                        bool reorder_space_filling)
{
  constexpr int DIM = 3;
  const auto range = bbox.range();

  auto* mesh =
    new mfem::Mesh(mfem::Mesh::MakeCartesian3D(res[0],
                                               res[1],
                                               res[2],
                                               mfem::Element::HEXAHEDRON,
                                               range[0],
                                               range[1],
                                               range[2],
                                               reorder_space_filling));

  // Offset the mesh to lie w/in the bounding box
  const int NV = mesh->GetNV();
  for(int i = 0; i < NV; ++i)
  {
    double* v = mesh->GetVertex(i);
    for(int d = 0; d < DIM; ++d)
    {
      v[d] += bbox.getMin()[d];
    }
  }

  // Ensure that mesh has high order nodes
  mesh->SetCurvature(polynomial_order);

  return mesh;
}

#endif  // AXOM_USE_MFEM

#if defined(AXOM_USE_SIDRE)

axom::sidre::Group* make_structured_blueprint_box_mesh(
  axom::sidre::Group* meshGrp,
  const primal::BoundingBox<double, 3>& bbox,
  const primal::NumericArray<int, 3>& res,
  const std::string& topologyName,
  const std::string& coordsetName,
  axom::runtime_policy::Policy runtimePolicy)
{
  auto* topoGrp = meshGrp->createGroup("topologies")->createGroup(topologyName);
  SLIC_ERROR_IF(topoGrp == nullptr,
                "Cannot allocate topology '" + topologyName +
                  "' in blueprint mesh '" + meshGrp->getName() +
                  "'.  It already exists.");

  auto* coordsetGrp =
    meshGrp->createGroup("coordsets")->createGroup(coordsetName);
  SLIC_ERROR_IF(coordsetGrp == nullptr,
                "Cannot allocate coordset '" + coordsetName +
                  "' in blueprint mesh '" + meshGrp->getName() +
                  "'.  It already exists.");

  topoGrp->createView("type")->setString("structured");
  topoGrp->createView("coordset")->setString(coordsetName);
  auto* dimsGrp = topoGrp->createGroup("elements/dims");

  constexpr int DIM = 3;

  auto ni = res[0];
  auto nj = res[1];
  auto nk = res[2];

  const axom::StackArray<axom::IndexType, DIM> vertsShape {res[0] + 1,
                                                           res[1] + 1,
                                                           res[2] + 1};
  const auto numVerts = vertsShape[0] * vertsShape[1] * vertsShape[2];

  dimsGrp->createViewScalar("i", ni);
  dimsGrp->createViewScalar("j", nj);
  dimsGrp->createViewScalar("k", nk);

  coordsetGrp->createView("type")->setString("explicit");
  auto* valuesGrp = coordsetGrp->createGroup("values");
  auto* xVu =
    valuesGrp->createViewAndAllocate("x",
                                     axom::sidre::DataTypeId::FLOAT64_ID,
                                     numVerts);
  auto* yVu =
    valuesGrp->createViewAndAllocate("y",
                                     axom::sidre::DataTypeId::FLOAT64_ID,
                                     numVerts);
  auto* zVu =
    valuesGrp->createViewAndAllocate("z",
                                     axom::sidre::DataTypeId::FLOAT64_ID,
                                     numVerts);

  const axom::MDMapping<DIM> vertMapping(vertsShape,
                                         axom::ArrayStrideOrder::COLUMN);
  axom::ArrayView<double, DIM> xView(xVu->getData(),
                                     vertsShape,
                                     vertMapping.strides());
  axom::ArrayView<double, DIM> yView(yVu->getData(),
                                     vertsShape,
                                     vertMapping.strides());
  axom::ArrayView<double, DIM> zView(zVu->getData(),
                                     vertsShape,
                                     vertMapping.strides());

  fill_cartesian_coords_3d(runtimePolicy, bbox, xView, yView, zView);

  #if defined(AXOM_DEBUG) && defined(AXOM_USE_CONDUIT)
  {
    conduit::Node info;
    bool isValid = verifyBlueprintMesh(meshGrp, info);
    SLIC_ASSERT_MSG(isValid, "Internal error: Generated mesh is invalid.");
  }
  #endif

  return meshGrp;
}

  #if defined(AXOM_USE_CONDUIT)
axom::sidre::Group* make_unstructured_blueprint_box_mesh(
  axom::sidre::Group* meshGrp,
  const primal::BoundingBox<double, 3>& bbox,
  const primal::NumericArray<int, 3>& res,
  const std::string& topologyName,
  const std::string& coordsetName,
  axom::runtime_policy::Policy runtimePolicy)
{
  make_structured_blueprint_box_mesh(meshGrp,
                                     bbox,
                                     res,
                                     topologyName,
                                     coordsetName,
                                     runtimePolicy);
  convert_blueprint_structured_explicit_to_unstructured(meshGrp, topologyName);
  return meshGrp;
}

void convert_blueprint_structured_explicit_to_unstructured(
  axom::sidre::Group* meshGrp,
  const std::string& topoName)
{
  const std::string& coordsetName =
    meshGrp->getView(axom::fmt::format("topologies/{}/coordset", topoName))
      ->getString();

  sidre::Group* coordsetGrp = nullptr;
    #if defined(AXOM_USE_UMPIRE)
  /* If using device memory, temporarily copy coords to host
     so we can use conduit's blueprint::mesh utilities.
     When we re-import the newCoords, we will get the data
     back to meshGrp memory space.
  */
  coordsetGrp = meshGrp->getGroup("coordsets")->getGroup(coordsetName);
  int coordsetAllocId = coordsetGrp->getDefaultAllocatorID();
  MemorySpace memSpace = detail::getAllocatorSpace(coordsetAllocId);
  sidre::Group* stashGrp = nullptr;
  sidre::Group* stashedValuesGrp = nullptr;
  if(memSpace == MemorySpace::Device)
  {
    int hostAllocId = execution_space<axom::SEQ_EXEC>::allocatorID();
    stashGrp = meshGrp->createGroup("tempStash");
    stashedValuesGrp = stashGrp->moveGroup(coordsetGrp->getGroup("values"));
    coordsetGrp->deepCopyGroup(stashedValuesGrp, hostAllocId);
  }
    #endif

  // Convert mesh to conduit::Node to use conduit's blueprint support.
  conduit::Node info;
  conduit::Node curMesh;
  meshGrp->createNativeLayout(curMesh);
  const conduit::Node& curTopo =
    curMesh.fetch_existing(axom::fmt::format("topologies/{}", topoName));
  SLIC_ASSERT(
    conduit::blueprint::mesh::topology::structured::verify(curTopo, info));

  // Use conduit to convert to unstructured.
  conduit::Node newTopo;
  conduit::Node newCoords;
  conduit::blueprint::mesh::topology::structured::to_unstructured(curTopo,
                                                                  newTopo,
                                                                  newCoords);

  // Copy unstructured back into meshGrp.
  meshGrp->getGroup("topologies")->destroyGroup(topoName);
  auto* topoGrp = meshGrp->getGroup("topologies")->createGroup(topoName);
  topoGrp->importConduitTree(newTopo);
  topoGrp->getView("coordset")
    ->setString(coordsetName);  // Is this needed?  Is coordset already set?

  meshGrp->getGroup("coordsets")->destroyGroup(coordsetName);
  coordsetGrp = meshGrp->getGroup("coordsets")->createGroup(coordsetName);
  coordsetGrp->importConduitTree(newCoords);

    #define ADD_EXTRA_DATA_FOR_MINT 1
    #if ADD_EXTRA_DATA_FOR_MINT
  /*
    Constructing a mint mesh from meshGrp fails unless we add some
    extra data.  Blueprint doesn't require this extra data.  (The mesh
    passes conduit's Blueprint verification.)  This should be fixed,
    or we should write better blueprint support utilities.
  */
  /*
    Make the coordinate arrays 2D to use mint::Mesh.
    For some reason, mint::Mesh requires the arrays to be
    2D, even though the second dimension is always 1.
  */
  axom::IndexType curShape[2];
  int curDim;
  auto* valuesGrp = coordsetGrp->getGroup("values");
  curDim = valuesGrp->getView("x")->getShape(2, curShape);
  assert(curDim == 1);
  axom::IndexType vertsShape[2] = {curShape[0], 1};
  valuesGrp->getView("x")->reshapeArray(2, vertsShape);
  valuesGrp->getView("y")->reshapeArray(2, vertsShape);
  valuesGrp->getView("z")->reshapeArray(2, vertsShape);

  // Make connectivity array 2D for the same reason.
  auto* elementsGrp = topoGrp->getGroup("elements");
  auto* connView = elementsGrp->getView("connectivity");
  curDim = connView->getShape(2, curShape);
  constexpr axom::IndexType NUM_VERTS_PER_HEX = 8;
  SLIC_ASSERT(curDim == 1);
  SLIC_ASSERT(curShape[0] % NUM_VERTS_PER_HEX == 0);
  axom::IndexType connShape[2] = {curShape[0] / NUM_VERTS_PER_HEX,
                                  NUM_VERTS_PER_HEX};
  connView->reshapeArray(2, connShape);

  // mint::Mesh requires connectivity strides, even though Blueprint doesn't.
  elementsGrp->createViewScalar("stride", NUM_VERTS_PER_HEX);

  // mint::Mesh requires field group, even though Blueprint doesn't.
  meshGrp->createGroup("fields");
    #endif

    #if defined(AXOM_DEBUG)
  conduit::Node newMesh;
  meshGrp->createNativeLayout(newMesh);
  SLIC_ASSERT(conduit::blueprint::mesh::verify(newMesh, info));
    #endif

  return;
}

bool verifyBlueprintMesh(const axom::sidre::Group* meshGrp, conduit::Node info)
{
  conduit::Node meshNode;
  meshGrp->createNativeLayout(meshNode);
  bool isValid = conduit::blueprint::mesh::verify(meshNode, info);
  if(!isValid) info.print();
  return isValid;
}
  #endif

#endif

void fill_cartesian_coords_3d(axom::runtime_policy::Policy runtimePolicy,
                              const primal::BoundingBox<double, 3>& domainBox,
                              axom::ArrayView<double, 3>& xView,
                              axom::ArrayView<double, 3>& yView,
                              axom::ArrayView<double, 3>& zView)
{
  if(runtimePolicy == axom::runtime_policy::Policy::seq)
  {
    fill_cartesian_coords_3d_impl<axom::SEQ_EXEC>(domainBox, xView, yView, zView);
  }
#if defined(AXOM_RUNTIME_POLICY_USE_OPENMP)
  if(runtimePolicy == axom::runtime_policy::Policy::omp)
  {
    fill_cartesian_coords_3d_impl<axom::OMP_EXEC>(domainBox, xView, yView, zView);
  }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_CUDA)
  if(runtimePolicy == axom::runtime_policy::Policy::cuda)
  {
    fill_cartesian_coords_3d_impl<axom::CUDA_EXEC<256>>(domainBox,
                                                        xView,
                                                        yView,
                                                        zView);
  }
#endif
#if defined(AXOM_RUNTIME_POLICY_USE_HIP)
  if(runtimePolicy == axom::runtime_policy::Policy::hip)
  {
    fill_cartesian_coords_3d_impl<axom::HIP_EXEC<256>>(domainBox,
                                                       xView,
                                                       yView,
                                                       zView);
  }
#endif
}

template <typename ExecSpace>
void fill_cartesian_coords_3d_impl(const primal::BoundingBox<double, 3>& domainBox,
                                   axom::ArrayView<double, 3>& xView,
                                   axom::ArrayView<double, 3>& yView,
                                   axom::ArrayView<double, 3>& zView)
{
  const auto& shape = xView.shape();
  const auto& mapping = xView.mapping();
  auto order = mapping.getStrideOrder();

  SLIC_ASSERT(shape == yView.shape());
  SLIC_ASSERT(shape == zView.shape());
  SLIC_ASSERT(mapping == yView.mapping());
  SLIC_ASSERT(mapping == zView.mapping());

  // Mesh resolution
  const axom::primal::NumericArray<int, 3> res {shape[0] - 1,
                                                shape[1] - 1,
                                                shape[2] - 1};

  // Mesh spacings.
  double dx = (domainBox.getMax()[0] - domainBox.getMin()[0]) / res[0];
  double dy = (domainBox.getMax()[1] - domainBox.getMin()[1]) / res[1];
  double dz = (domainBox.getMax()[2] - domainBox.getMin()[2]) / res[2];

#if defined(AXOM_USE_RAJA)
  RAJA::RangeSegment kRange(0, shape[2]);
  RAJA::RangeSegment jRange(0, shape[1]);
  RAJA::RangeSegment iRange(0, shape[0]);
  using EXEC_POL =
    typename axom::internal::nested_for_exec<ExecSpace>::loop3d_policy;
  if(int(order) & int(axom::ArrayStrideOrder::COLUMN))
  {
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(iRange, jRange, kRange),
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j, axom::IndexType k) {
        xView(i, j, k) = domainBox.getMin()[0] + i * dx;
        yView(i, j, k) = domainBox.getMin()[1] + j * dy;
        zView(i, j, k) = domainBox.getMin()[2] + k * dz;
      });
  }
  else
  {
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(kRange, jRange, iRange),
      AXOM_LAMBDA(axom::IndexType k, axom::IndexType j, axom::IndexType i) {
        xView(i, j, k) = domainBox.getMin()[0] + i * dx;
        yView(i, j, k) = domainBox.getMin()[1] + j * dy;
        zView(i, j, k) = domainBox.getMin()[2] + k * dz;
      });
  }
#else
  for(int k = 0; k < res[2]; ++k)
  {
    for(int j = 0; j < res[1]; ++j)
    {
      for(int i = 0; i < res[0]; ++i)
      {
        xView(i, j, k) = domainBox.getMin()[0] + i * dx;
        yView(i, j, k) = domainBox.getMin()[1] + j * dy;
        zView(i, j, k) = domainBox.getMin()[2] + k * dz;
      }
    }
  }
#endif

  return;
}

}  // namespace util
}  // namespace quest
}  // namespace axom
