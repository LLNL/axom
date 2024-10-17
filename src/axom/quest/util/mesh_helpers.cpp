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
#include "axom/core/WhereMacro.hpp"

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
  const std::string& coordsetName)
{
  constexpr int DIM = 3;

  auto ni = res[0];
  auto nj = res[1];
  auto nk = res[2];

  const axom::StackArray<axom::IndexType, DIM> vertsShape {res[0] + 1,
                                                           res[1] + 1,
                                                           res[2] + 1};
  const auto numVerts = vertsShape[0] * vertsShape[1] * vertsShape[2];

  auto* topoGrp = meshGrp->createGroup("topologies")->createGroup(topologyName);
  auto* coordsetGrp =
    meshGrp->createGroup("coordsets")->createGroup(coordsetName);

  topoGrp->createView("type")->setString("structured");
  topoGrp->createView("coordset")->setString(coordsetName);
  auto* dimsGrp = topoGrp->createGroup("elements/dims");
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

  const axom::MDMapping<DIM> vertMapping(vertsShape, axom::ArrayStrideOrder::ROW);
  axom::ArrayView<double, DIM> xView(xVu->getData(),
                                     vertsShape,
                                     vertMapping.strides());
  axom::ArrayView<double, DIM> yView(yVu->getData(),
                                     vertsShape,
                                     vertMapping.strides());
  axom::ArrayView<double, DIM> zView(zVu->getData(),
                                     vertsShape,
                                     vertMapping.strides());

  double dx = bbox.getMax()[0] - bbox.getMin()[0];
  double dy = bbox.getMax()[1] - bbox.getMin()[1];
  double dz = bbox.getMax()[2] - bbox.getMin()[2];
  for(int k = 0; k < nk + 1; ++k)
  {
    for(int j = 0; j < nj + 1; ++j)
    {
      for(int i = 0; i < ni + 1; ++i)
      {
        xView(i, j, k) = bbox.getMin()[0] + i * dx;
        yView(i, j, k) = bbox.getMin()[1] + j * dy;
        zView(i, j, k) = bbox.getMin()[2] + k * dz;
      }
    }
  }

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
  const std::string& coordsetName)
{
  make_structured_blueprint_box_mesh(meshGrp, bbox, res, topologyName, coordsetName);
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
  auto* coordsetGrp = meshGrp->getGroup("coordsets")->createGroup(coordsetName);
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

}  // namespace util
}  // namespace quest
}  // namespace axom
