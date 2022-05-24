// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file scattered_interpolation.cpp
 * \brief This examples uses a Delaunay mesh to perform scattered data
 * interpolation on fields defined over a point set
 */

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/sidre.hpp"
#include "axom/primal.hpp"
#include "axom/quest.hpp"

#include "conduit_blueprint.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <memory>

namespace primal = axom::primal;
namespace quest = axom::quest;
namespace sidre = axom::sidre;

namespace internal
{
namespace blueprint
{
/**
 *  \brief Simple wrapper to a blueprint particle mesh
 *
 *  Given a sidre Group, creates the stubs for a mesh blueptint particle mesh
 */
struct PointMesh
{
public:
  explicit PointMesh(sidre::Group* group = nullptr,
                     const std::string& coordset = "coords",
                     const std::string& topology = "mesh")
    : m_group(group)
    , m_meshName(topology)
    , m_coordName(coordset)
  {
    setBlueprintGroup(m_group, coordset, topology);
  }
  /// Gets the root group for this mesh blueprint
  sidre::Group* rootGroup() const { return m_group; }
  /// Gets the parent group for the blueprint coordinate set
  sidre::Group* coordsGroup() const { return m_coordsGroup; }
  /// Gets the parent group for the blueprint mesh topology
  sidre::Group* topoGroup() const { return m_topoGroup; }

  /// Returns true if points have been added to the particle mesh
  bool hasPoints() const
  {
    return m_coordsGroup != nullptr && m_coordsGroup->hasGroup("values");
  }

  /// Returns the number of points in the particle mesh
  int numPoints() const
  {
    return hasPoints() ? m_coordsGroup->getView("values/x")->getNumElements() : 0;
  }

  int dimension() const
  {
    if(hasPoints())
    {
      return m_coordsGroup->hasView("values/z")
        ? 3
        : (m_coordsGroup->hasView("values/y") ? 2 : 1);
    }
    return 0;
  }

  /**
   * Sets the parent group for the entire mesh and sets up the blueprint stubs
   * for the "coordset", "topologies", "fields" and "state"
   */
  void setBlueprintGroup(sidre::Group* group,
                         const std::string& coordset = "coords",
                         const std::string& topology = "mesh")
  {
    // TODO: Ensure that we delete previous hierarchy if it existed

    m_group = group;
    m_meshName = topology;
    m_coordName = coordset;

    if(m_group != nullptr)
    {
      createBlueprintStubs(coordset, topology);
    }
  }

  /// Set the coordinate data from an array of primal Points, templated on the dimension
  template <int NDIMS>
  void setPoints(const axom::Array<primal::Point<double, NDIMS>>& pts)
  {
    SLIC_ASSERT_MSG(m_group != nullptr,
                    "Must set blueprint group before setPoints()");

    const int SZ = pts.size();

    // create views into a shared buffer for the coordinates, with stride NDIMS
    auto* buf =
      m_group->getDataStore()->createBuffer(sidre::DOUBLE_ID, NDIMS * SZ)->allocate();
    switch(NDIMS)
    {
    case 3:
      m_coordsGroup->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, NDIMS);
      m_coordsGroup->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, NDIMS);
      m_coordsGroup->createView("values/z")->attachBuffer(buf)->apply(SZ, 2, NDIMS);
      break;
    case 2:
      m_coordsGroup->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, NDIMS);
      m_coordsGroup->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, NDIMS);
      break;
    default:
      m_coordsGroup->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, NDIMS);
      break;
    }

    // copy coordinate data into the buffer
    const std::size_t nbytes = sizeof(double) * SZ * NDIMS;
    axom::copy(buf->getVoidPtr(), pts.data(), nbytes);

    // set the default connectivity
    sidre::Array<int> arr(m_topoGroup->createView("elements/connectivity"), SZ, SZ);
    for(int i = 0; i < SZ; ++i)
    {
      arr[i] = i;
    }
  }

  template <typename T>
  void registerNodalScalarField(const std::string& fieldName)
  {
    SLIC_ASSERT_MSG(hasPoints(),
                    "Cannot register a field with the BlueprintParticleMesh "
                    "before adding points");

    auto* fld = m_fieldsGroup->createGroup(fieldName);
    fld->createViewString("association", "vertex");
    fld->createViewString("topology", m_topoGroup->getName());
    fld->createViewAndAllocate("values",
                               sidre::detail::SidreTT<T>::id,
                               numPoints());
  }

  template <typename T>
  void registerNodalVectorField(const std::string& fieldName)
  {
    SLIC_ASSERT_MSG(hasPoints(),
                    "Cannot register a field with the BlueprintParticleMesh "
                    "before adding points");

    const int SZ = numPoints();
    const int DIM = dimension();

    auto* fld = m_fieldsGroup->createGroup(fieldName);
    fld->createViewString("association", "vertex");
    fld->createViewString("topology", m_topoGroup->getName());

    // create views into a shared buffer for the coordinates, with stride NDIMS
    auto* buf = m_group->getDataStore()
                  ->createBuffer(sidre::detail::SidreTT<T>::id, DIM * SZ)
                  ->allocate();
    switch(DIM)
    {
    case 3:
      fld->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, DIM);
      fld->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, DIM);
      fld->createView("values/z")->attachBuffer(buf)->apply(SZ, 2, DIM);
      break;
    case 2:
      fld->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, DIM);
      fld->createView("values/y")->attachBuffer(buf)->apply(SZ, 1, DIM);
      break;
    default:
      fld->createView("values/x")->attachBuffer(buf)->apply(SZ, 0, DIM);
      break;
    }
  }

  bool hasField(const std::string& fieldName) const
  {
    return m_fieldsGroup->hasGroup(fieldName);
  }

  template <typename T>
  axom::ArrayView<T> getNodalScalarField(const std::string& fieldName) const
  {
    SLIC_ASSERT_MSG(hasPoints(),
                    "Cannot extract a field from the BlueprintParticleMesh "
                    "before adding points");

    T* data = hasField(fieldName)
      ? static_cast<T*>(
          m_fieldsGroup->getView(axom::fmt::format("{}/values", fieldName))
            ->getVoidPtr())
      : nullptr;

    return axom::ArrayView<T>(data, numPoints());
  }

  template <typename T>
  axom::ArrayView<T> getNodalVectorField(const std::string& fieldName) const
  {
    SLIC_ASSERT_MSG(hasPoints(),
                    "Cannot extract a field from the BlueprintParticleMesh "
                    "before adding points");

    T* data = hasField(fieldName)
      ? static_cast<T*>(
          m_fieldsGroup->getView(axom::fmt::format("{}/values/x", fieldName))
            ->getVoidPtr())
      : nullptr;

    return axom::ArrayView<T>(data, numPoints());
  }

  /// Checks whether the blueprint is valid and prints diagnostics
  bool isValid() const
  {
    conduit::Node mesh_node;
    m_group->createNativeLayout(mesh_node);

    bool success = true;
    conduit::Node info;
    if(!conduit::blueprint::verify("mesh", mesh_node, info))
    {
      SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
      success = false;
    }

    return success;
  }

  /// Outputs the mesh to disk
  void saveMesh(const std::string& outputMesh, const std::string& protocol)
  {
    auto* ds = m_group->getDataStore();

    const std::string ext = (protocol.find("hdf5") != std::string::npos)
      ? "hdf5"
      : (protocol.find("json") != std::string::npos ? "json" : protocol);
    const std::string fname = axom::fmt::format("{}.{}", outputMesh, ext);

    // generate metadata and write out as "root" file
    {
      auto* rootfile_grp = ds->getRoot()->hasGroup("rootfile_data")
        ? ds->getRoot()->getGroup("rootfile_data")
        : ds->getRoot()->createGroup("rootfile_data");

      rootfile_grp->createViewScalar("number_of_files", 1);
      rootfile_grp->createViewString("file_pattern", fname);
      rootfile_grp->createViewScalar("number_of_trees", 1);
      rootfile_grp->createViewString("tree_pattern", "/");
      rootfile_grp->createViewString("protocol/name", protocol);
      rootfile_grp->createViewString("protocol/version", "0.8.2");

      ds->generateBlueprintIndex(
        m_group->getPathName(),
        m_group->getName(),
        axom::fmt::format("rootfile_data/blueprint_index/{}", m_group->getName()),
        1);

      rootfile_grp->save(axom::fmt::format("{}.root", outputMesh), "json");

      ds->getRoot()->destroyGroup("rootfile_data");
    }

    ds->getRoot()->save(fname, protocol);
  }

private:
  /// Creates blueprint stubs for this mesh
  void createBlueprintStubs(const std::string& coords, const std::string& topo)
  {
    SLIC_ASSERT(m_group != nullptr);

    m_coordsGroup = m_group->createGroup("coordsets")->createGroup(coords);
    m_coordsGroup->createViewString("type", "explicit");
    m_coordsGroup->createGroup("values");

    m_topoGroup = m_group->createGroup("topologies")->createGroup(topo);
    m_topoGroup->createViewString("coordset", coords);
    m_topoGroup->createViewString("type", "unstructured");
    m_topoGroup->createViewString("elements/shape", "point");

    m_fieldsGroup = m_group->createGroup("fields");
  }

private:
  sidre::Group* m_group;

  sidre::Group* m_coordsGroup;
  sidre::Group* m_topoGroup;
  sidre::Group* m_fieldsGroup;

  std::string m_meshName;
  std::string m_coordName;
};

}  // namespace blueprint
}  // namespace internal

/// Struct to parse and contain command line arguments
struct Input
{
  std::string outputFile {"scattered_interpolation"};
  int numRandPoints {20};
  int numQueryPoints {20};
  int dimension {2};
  std::vector<double> boundsMin;
  std::vector<double> boundsMax;
  std::string outputProtocol = SIDRE_DEFAULT_PROTOCOL;

  const std::set<std::string> s_validProtocols {"json",
                                                "sidre_json"
#ifdef AXOM_USE_HDF5
                                                ,
                                                "sidre_hdf5"
#endif
  };

public:
  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-n,--nrandpt", numRandPoints)
      ->description("The number of points in the input mesh")
      ->capture_default_str();

    app.add_option("-q,--nquerypt", numQueryPoints)
      ->description("The number of query points")
      ->capture_default_str();

    app.add_option("-d,--dim", dimension)
      ->description(
        "The dimension of the mesh. 2 for triangle mesh in 2D; "
        "3 for tetrahedral mesh in 3D")
      ->capture_default_str();

    app.add_option("-o,--outfile", outputFile)
      ->description("The output file")
      ->capture_default_str();

    app.add_option("-p,--protocol", outputProtocol)
      ->description("Set the output protocol for sidre point meshes")
      ->capture_default_str()
      ->check(axom::CLI::IsMember(s_validProtocols));

    // Optional bounding box for query region
    auto* minbb = app.add_option("--min", boundsMin)
                    ->description("Min bounds for query box (x,y[,z])")
                    ->expected(2, 3);
    auto* maxbb = app.add_option("--max", boundsMax)
                    ->description(
                      "Max bounds for query box (x,y[,z]). "
                      "Defaults to unit square/cube when not provided")
                    ->expected(2, 3);
    minbb->needs(maxbb);
    maxbb->needs(minbb);

    app.get_formatter()->column_width(45);

    // could throw an exception
    app.parse(argc, argv);

    // If user doesn't provide bounds, default to unit cube
    if(boundsMin.empty())
    {
      boundsMin.clear();
      boundsMax.clear();
      for(int i = 0; i < dimension; ++i)
      {
        boundsMin.push_back(0.);
        boundsMax.push_back(1.);
      }
    }

    SLIC_INFO(axom::fmt::format(R"(Using parameter values:
    {{
      dimension: {}
      nrandpt: {}
      nquerypt: {}
      bounding box min: {{{}}}
      bounding box max: {{{}}}
      outfile = '{}'
      output protocol = '{}'
    }})",
                                dimension,
                                numRandPoints,
                                numQueryPoints,
                                axom::fmt::join(boundsMin, ", "),
                                axom::fmt::join(boundsMax, ", "),
                                outputFile,
                                outputProtocol));
  }
};

template <int DIM>
axom::Array<primal::Point<double, DIM>> generatePts(
  int numPts,
  const std::vector<double>& bb_min,
  const std::vector<double>& bb_max)
{
  using axom::utilities::random_real;

  using PointType = typename primal::Point<double, DIM>;
  using BoundingBox = typename primal::BoundingBox<double, DIM>;

  axom::Array<PointType> pts(numPts, numPts);

  BoundingBox bbox {PointType(bb_min.data()), PointType(bb_max.data())};

  // generate random points within bounding box
  for(int i = 0; i < numPts; ++i)
  {
    PointType& pt = pts[i];
    for(int d = 0; d < DIM; ++d)
    {
      pt[d] = random_real(bbox.getMin()[d], bbox.getMax()[d]);
    }
  }

  return pts;
}

int main(int argc, char** argv)
{
  // Initialize the SLIC logger
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  // Initialize default parameters and update with command line arguments:
  Input params;
  axom::CLI::App app {"Scattered interpolation over 2D or 3D point collection"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

  sidre::DataStore ds;

  const std::string input_coords_name = "input_coords";
  const std::string input_mesh_name = "input_mesh";
  auto* input_pt_mesh_group = ds.getRoot()->createGroup("input_point_mesh");
  internal::blueprint::PointMesh inputMesh(input_pt_mesh_group,
                                           input_coords_name,
                                           input_mesh_name);

  const std::string query_coords_name = "query_coords";
  const std::string query_mesh_name = "query_mesh";
  auto* query_pt_mesh_group = ds.getRoot()->createGroup("query_point_mesh");
  internal::blueprint::PointMesh queryMesh(query_pt_mesh_group,
                                           query_coords_name,
                                           query_mesh_name);

  // Initialize the mesh points
  switch(params.dimension)
  {
  case 2:
    inputMesh.setPoints(
      generatePts<2>(params.numRandPoints, params.boundsMin, params.boundsMax));

    queryMesh.setPoints(
      generatePts<2>(params.numQueryPoints, params.boundsMin, params.boundsMax));
    break;
  case 3:
    inputMesh.setPoints(
      generatePts<3>(params.numRandPoints, params.boundsMin, params.boundsMax));

    queryMesh.setPoints(
      generatePts<3>(params.numQueryPoints, params.boundsMin, params.boundsMax));
    break;
  }

  // Add a position field for our interpolation
  switch(params.dimension)
  {
  case 2:
  {
    using PtType = primal::Point<double, 2>;
    inputMesh.registerNodalScalarField<double>("pos_x");
    inputMesh.registerNodalScalarField<double>("pos_y");
    auto pos_x = inputMesh.getNodalScalarField<double>("pos_x");
    auto pos_y = inputMesh.getNodalScalarField<double>("pos_y");

    const int nPts = inputMesh.numPoints();
    auto pos = axom::ArrayView<PtType>(
      static_cast<PtType*>(
        inputMesh.coordsGroup()->getView("values/x")->getVoidPtr()),
      nPts);

    for(int i = 0; i < nPts; ++i)
    {
      auto& pt = pos[i];
      pos_x[i] = pt[0];
      pos_y[i] = pt[1];
    }
  }
  break;
  case 3:
  {
    using PtType = primal::Point<double, 3>;
    inputMesh.registerNodalScalarField<double>("pos_x");
    inputMesh.registerNodalScalarField<double>("pos_y");
    inputMesh.registerNodalScalarField<double>("pos_z");
    auto pos_x = inputMesh.getNodalScalarField<double>("pos_x");
    auto pos_y = inputMesh.getNodalScalarField<double>("pos_y");
    auto pos_z = inputMesh.getNodalScalarField<double>("pos_z");

    const int nPts = inputMesh.numPoints();
    auto pos = axom::ArrayView<PtType>(
      static_cast<PtType*>(
        inputMesh.coordsGroup()->getView("values/x")->getVoidPtr()),
      nPts);

    for(int i = 0; i < nPts; ++i)
    {
      auto& pt = pos[i];
      pos_x[i] = pt[0];
      pos_y[i] = pt[1];
      pos_z[i] = pt[2];
    }
  }
  break;
  }

  queryMesh.registerNodalScalarField<axom::IndexType>("cell_idx");
  queryMesh.registerNodalScalarField<double>("pos_x");
  queryMesh.registerNodalScalarField<double>("pos_y");
  if(params.dimension > 2)
  {
    queryMesh.registerNodalScalarField<double>("pos_z");
  }

  // Write input mesh to file
  {
    std::string file = params.outputFile + "_mesh_after_generate";
    SLIC_INFO(
      axom::fmt::format("After generate points -- writing mesh file: '{}'", file));
    inputMesh.saveMesh(file, params.outputProtocol);
  }

  std::unique_ptr<quest::ScatteredInterpolation<2>> scattered_2d;
  std::unique_ptr<quest::ScatteredInterpolation<3>> scattered_3d;

  // Convert blueprint mesh from Sidre to Conduit
  conduit::Node bp_input;
  inputMesh.rootGroup()->createNativeLayout(bp_input);

  conduit::Node bp_query;
  queryMesh.rootGroup()->createNativeLayout(bp_query);

  // Initialize the mesh points
  switch(params.dimension)
  {
  case 2:
    scattered_2d = std::unique_ptr<quest::ScatteredInterpolation<2>>(
      new quest::ScatteredInterpolation<2>);
    scattered_2d->buildTriangulation(bp_input, input_coords_name);
    scattered_2d->locatePoints(bp_query, query_coords_name);

    scattered_2d->interpolateField(bp_query,
                                   query_coords_name,
                                   bp_input,
                                   "pos_x",
                                   "pos_x");
    scattered_2d->interpolateField(bp_query,
                                   query_coords_name,
                                   bp_input,
                                   "pos_y",
                                   "pos_y");
    break;
  case 3:
    scattered_3d = std::unique_ptr<quest::ScatteredInterpolation<3>>(
      new quest::ScatteredInterpolation<3>);
    scattered_3d->buildTriangulation(bp_input, input_coords_name);
    scattered_3d->locatePoints(bp_query, query_coords_name);

    scattered_3d->interpolateField(bp_query,
                                   query_coords_name,
                                   bp_input,
                                   "pos_x",
                                   "pos_x");
    scattered_3d->interpolateField(bp_query,
                                   query_coords_name,
                                   bp_input,
                                   "pos_y",
                                   "pos_y");
    scattered_3d->interpolateField(bp_query,
                                   query_coords_name,
                                   bp_input,
                                   "pos_z",
                                   "pos_z");
    break;
  }

  // Write query mesh to file
  {
    std::string file = params.outputFile + "_mesh_after_query";
    SLIC_INFO(
      axom::fmt::format("After locating points on Delaunay triangulation: '{}'",
                        file));
    queryMesh.saveMesh(file, params.outputProtocol);
  }

  return 0;
}
