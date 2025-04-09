// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/mint.hpp"
#include "axom/slam.hpp"
#include "axom/spin.hpp"
#include "axom/quest.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include <string>

namespace slam = axom::slam;
namespace primal = axom::primal;
namespace quest = axom::quest;
namespace spin = axom::spin;

template <int TOPO_DIM, int SPACE_DIM = 3>
struct SimplicialMesh
{
public:
  static constexpr int NODES_PER_ZONE = TOPO_DIM + 1;

  using PositionType = slam::DefaultPositionType;

  // types for sets
  using ElementType = slam::DefaultElementType;
  using PointType = primal::Point<double, SPACE_DIM>;
  using NodeSet = slam::PositionSet<PositionType, ElementType>;
  using ZoneSet = slam::PositionSet<PositionType, ElementType>;

  // types for relations
  using Indir = slam::policies::ArrayIndirection<PositionType, ElementType>;
  using ZNStride = slam::policies::CompileTimeStride<PositionType, NODES_PER_ZONE>;
  using ConstantCardinality = slam::policies::ConstantCardinality<PositionType, ZNStride>;
  using ZoneToNodeRelation =
    typename slam::StaticRelation<PositionType, ElementType, ConstantCardinality, Indir, ZoneSet, NodeSet>;
  using ZNBuilder = typename ZoneToNodeRelation::RelationBuilder;

  /// types for maps
  using Vertices =
    slam::Map<PointType, NodeSet, slam::policies::ArrayIndirection<PositionType, PointType>>;

public:
  SimplicialMesh(int num_verts = 0, int num_zones = 0) { reset(num_verts, num_zones); }

  void reset(int num_verts = 0, int num_zones = 0)
  {
    nodeSet = NodeSet(num_verts);
    zoneSet = ZoneSet(num_zones);

    // set up the zone-node relation
    m_zn.resize(zoneSet.size() * NODES_PER_ZONE);
    znRel = ZNBuilder().fromSet(&zoneSet).toSet(&nodeSet).indices(
      typename ZNBuilder::IndicesSetBuilder().size(m_zn.size()).data(&m_zn));

    // set up the vertex coords
    coords = Vertices(&nodeSet);

    this->checkValid();
  }

  bool checkValid(bool verboseOutput = false) const
  {
    SLIC_ASSERT(nodeSet.isValid(verboseOutput));
    SLIC_ASSERT(zoneSet.isValid(verboseOutput));
    if(nodeSet.size() > 0 && zoneSet.size() > 0)
    {
      SLIC_ASSERT(znRel.isValid(verboseOutput));
    }
    SLIC_ASSERT(coords.isValid(verboseOutput));

    return true;
  }

  void reindexMesh(const axom::Array<PositionType>& remap_verts, int uniqueVertCount)
  {
    // update node set and vertices
    {
      NodeSet newNodeSet(uniqueVertCount);
      Vertices newCoords(&newNodeSet);
      PositionType curr = 0;
      for(auto n : nodeSet)
      {
        if(remap_verts[n] == n)
        {
          newCoords[curr++] = coords[n];
        }
      }
      nodeSet = newNodeSet;
      coords = newCoords;
    }

    // no change to zoneSet

    // update zn data (we can do it directly rather than go through the relation)
    for(int i = 0; i < m_zn.size(); ++i)
    {
      m_zn[i] = remap_verts[m_zn[i]];
    }

    // update zn relation
    znRel = ZNBuilder().fromSet(&zoneSet).toSet(&nodeSet).indices(
      typename ZNBuilder::IndicesSetBuilder().size(m_zn.size()).data(&m_zn));

    checkValid();
  }

public:
  NodeSet nodeSet;
  ZoneSet zoneSet;

  ZoneToNodeRelation znRel;
  Vertices coords;

  axom::Array<PositionType> m_zn;
};

using TetMesh = SimplicialMesh<3, 3>;  // Volumetric tetrahedral mesh in 3D
using TriMesh = SimplicialMesh<2, 3>;  // Surface triangle mesh in 3D

/// Loads the pro-e tet mesh and returns a TetMesh instance
/// Fixes orientations of the tets, but defers fixing degeneracies until the vertices have been welded
TetMesh loadProe(const std::string& filename, bool verbose)
{
  namespace mint = axom::mint;

  SLIC_INFO(axom::fmt::format("Reading pro-e file: '{}' ...", filename));

  /// read pro-e mesh into a mint unstructured mesh
  auto mint_mesh = std::make_unique<mint::UnstructuredMesh<mint::SINGLE_SHAPE>>(3, axom::mint::TET);
  {
    axom::quest::ProEReader reader;
    reader.setFileName(filename);

    reader.read();
    SLIC_INFO(axom::fmt::format("Input mesh has {} vertices and {} tets",
                                reader.getNumNodes(),
                                reader.getNumTets()));

    reader.getMesh(mint_mesh.get());
  }

  if(mint_mesh == nullptr)
  {
    SLIC_ERROR("Mesh is null after reading the proe file.");
    return TetMesh();
  }

  if(verbose)
  {
    std::string vtk_output_file = "proe_tet.vtk";
    mint::write_vtk(mint_mesh.get(), vtk_output_file);
    SLIC_INFO("Saved tet mesh to VTK file: " << vtk_output_file);
  }

  /// convert from mint to our slam-based TetMesh
  TetMesh slam_mesh(mint_mesh->getNumberOfNodes(), mint_mesh->getNumberOfCells());

  // copy vertex data
  {
    const double* x = mint_mesh->getCoordinateArray(mint::X_COORDINATE);
    const double* y = mint_mesh->getCoordinateArray(mint::Y_COORDINATE);
    const double* z = mint_mesh->getCoordinateArray(mint::Z_COORDINATE);
    for(auto n : slam_mesh.nodeSet)
    {
      slam_mesh.coords[n] = TetMesh::PointType {x[n], y[n], z[n]};
    }
  }

  // copy tet-vertex connectivity data
  {
    using SpatialTet = primal::Tetrahedron<double, 3>;
    using IndexTet = decltype(slam_mesh.znRel[0]);

    auto spatialTet = [&](IndexTet index_tet) -> SpatialTet {
      return SpatialTet(slam_mesh.coords[index_tet[0]],
                        slam_mesh.coords[index_tet[1]],
                        slam_mesh.coords[index_tet[2]],
                        slam_mesh.coords[index_tet[3]]);
    };

    for(auto t : slam_mesh.zoneSet)
    {
      const axom::IndexType* connec = mint_mesh->getCellNodeIDs(t);
      auto tet = slam_mesh.znRel[t];
      tet[0] = connec[0];
      tet[1] = connec[1];
      tet[2] = connec[2];
      tet[3] = connec[3];

      // fix orientations -- tets should not have negative volumes
      if(spatialTet(tet).signedVolume() < 0.)
      {
        axom::utilities::swap(tet[2], tet[3]);
      }
    }
  }

  slam_mesh.checkValid(verbose);

  return slam_mesh;
}

/// Welds all vertices in the tet mesh that are within \a weldThresh of each other
void weldVertices(TetMesh& mesh, double weldThresh, bool verbose)
{
  SLIC_ASSERT_MSG(weldThresh > 0.,
                  "Welding threshold must be greater than 0. Passed in value was " << weldThresh);

  /// compute bounding box of mesh
  axom::primal::BoundingBox<double, 3> meshBB;
  for(auto n : mesh.nodeSet)
  {
    meshBB.addPoint(mesh.coords[n]);
  }
  SLIC_INFO(axom::fmt::format("Bounding box of mesh: {}", meshBB));

  /// create spatial index
  // Note: Use 64-bit index to accomodate small values of epsilon
  using IdxType = std::int64_t;
  using Lattice3 = spin::RectangularLattice<3, double, IdxType>;
  using GridCell = Lattice3::GridCell;

  // Define a lambda for hashing points (implementation from boost's hash_combine())
  auto point_hash = [](const GridCell& pt) {
    auto seed = std::hash<IdxType> {}(pt[0]);
    for(int i = 1; i < GridCell::DIMENSION; ++i)
    {
      seed ^= std::hash<IdxType> {}(pt[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  };

  // Run welding algorithm twice -- on the original grid and on a slightly translated grid
  // This ensures that the grid resolution doesn't affect the results
  for(const auto& offset : {TetMesh::PointType(0.), TetMesh::PointType(weldThresh / 2.)})
  {
    // We're going to find the unique indices for the welded vertices
    const int numVerts = mesh.nodeSet.size();
    int uniqueVertCount = 0;

    // Set up the lattice for quantizing points to an integer lattice
    TetMesh::PointType origin(meshBB.getMin().array() - offset.array());
    Lattice3 lattice(origin, Lattice3::SpaceVector(TetMesh::PointType(weldThresh)));

    // A hashmap from GridCells to the new vertex indices
    using GridCellToIndexMap = std::unordered_map<GridCell, IdxType, decltype(point_hash)>;
    GridCellToIndexMap vertexIndexMap(numVerts, point_hash);

    axom::Array<TetMesh::PositionType> vertex_remap(numVerts);
    for(auto n : mesh.nodeSet)
    {
      // find the new vertex index; if not already present, we'll keep this vertex
      auto res =
        vertexIndexMap.insert(std::make_pair(lattice.gridCell(mesh.coords[n]), uniqueVertCount));
      if(res.second == true)
      {
        uniqueVertCount++;
      }
      vertex_remap[n] = static_cast<TetMesh::PositionType>(res.first->second);
    }

    // Finally, update the mesh using the new vertex mapping
    mesh.reindexMesh(vertex_remap, uniqueVertCount);

    SLIC_INFO_IF(
      verbose,
      axom::fmt::format("After reindexing with offset = {}, mesh has {} vertices and {} tets",
                        offset[0],
                        mesh.nodeSet.size(),
                        mesh.zoneSet.size()));
  }
}

/// Creates a \a tri_mesh out of the boundary faces of the \a tet_mesh
void extractBoundaryFaces(const TetMesh& tet_mesh, TriMesh& tri_mesh, bool verbose)
{
  using PositionType = TetMesh::PositionType;
  using IndexTriple = std::tuple<PositionType, PositionType, PositionType>;
  using ZoneFacePair = std::pair<PositionType, int>;
  std::map<IndexTriple, ZoneFacePair> unmatched_faces;

  auto is_degenerate = [](const auto& tet) -> bool {
    return tet[0] == tet[1] || tet[0] == tet[2] || tet[0] == tet[3]  //
      || tet[1] == tet[2] || tet[1] == tet[3]                        //
      || tet[2] == tet[3];
  };

  // we sort the vertex indices to avoid order dependence,
  // but keep the association to the original tet face so we don't mess up the orientations
  auto create_tuple = [](PositionType idx0, PositionType idx1, PositionType idx2, int face) {
    axom::utilities::detail::ifswap(idx0, idx1);
    axom::utilities::detail::ifswap(idx1, idx2);
    axom::utilities::detail::ifswap(idx0, idx1);
    return std::make_pair(std::make_tuple(idx0, idx1, idx2), face);
  };

  for(auto t : tet_mesh.zoneSet)
  {
    auto tet = tet_mesh.znRel[t];

    if(is_degenerate(tet))
    {
      continue;
    }

    // assume input tet mesh is a manifold w/ boundaries
    // internal faces will appear twice; boundary faces once
    for(auto triple_face : {create_tuple(tet[1], tet[2], tet[3], 0),
                            create_tuple(tet[0], tet[2], tet[3], 1),
                            create_tuple(tet[0], tet[1], tet[3], 2),
                            create_tuple(tet[0], tet[1], tet[2], 3)})
    {
      const auto it = unmatched_faces.find(triple_face.first);
      if(it != unmatched_faces.end())
      {
        // SLIC_INFO(axom::fmt::format("Removing triple {}", triple_face.first));
        unmatched_faces.erase(it);
      }
      else
      {
        // SLIC_INFO(axom::fmt::format("Adding triple {}", triple_face.first));
        unmatched_faces.insert({triple_face.first, std::make_pair(t, triple_face.second)});
      }
    }
  }
  SLIC_INFO(axom::fmt::format("Retained {} boundary faces, out of {}",
                              unmatched_faces.size(),
                              tet_mesh.znRel.totalSize()));

  /// Copy the boundary faces into the triangle mesh
  const int num_verts = tet_mesh.nodeSet.size();
  const int num_triangles = unmatched_faces.size();
  tri_mesh.reset(num_verts, num_triangles);

  SLIC_INFO_IF(verbose,
               axom::fmt::format(
                 "Triangle mesh has {} vertices and {} triangles, for a total of {} tv entities",
                 tri_mesh.nodeSet.size(),
                 tri_mesh.zoneSet.size(),
                 tri_mesh.znRel.totalSize()));

  // 1. copy the vertices
  // Since STL doesn't care about vertex indices, our tri mesh will have all vertices from the tet mesh
  // Alternatively, we could process this a bit more and keep only the referenced boundary vertices
  for(auto n : tet_mesh.nodeSet)
  {
    tri_mesh.coords[n] = tet_mesh.coords[n];
  }

  // 2. copy the boundary triangles
  int curr = 0;
  for(const auto& kv : unmatched_faces)
  {
    // SLIC_INFO_IF(verbose, axom::fmt::format("Adding face {} for triple {}", curr, kv.first));

    const auto zone = kv.second.first;
    auto tet = tet_mesh.znRel[zone];

    const int face = kv.second.second;
    auto tri = tri_mesh.znRel[curr++];

    // Note -- the triple was sorted, so we go back to the indices from the face's owning tet
    switch(face)
    {
    case 0:
      tri[0] = tet[1];
      tri[1] = tet[2];
      tri[2] = tet[3];
      break;
    case 1:
      tri[0] = tet[0];
      tri[1] = tet[3];
      tri[2] = tet[2];
      break;
    case 2:
      tri[0] = tet[0];
      tri[1] = tet[1];
      tri[2] = tet[3];
      break;
    case 3:
      tri[0] = tet[0];
      tri[1] = tet[2];
      tri[2] = tet[1];
      break;
    }
  }
}

/// Write the triangle mesh out as an ascii STL mesh
bool writeSTL(const TriMesh& tri_mesh, const std::string& filename)
{
  axom::fmt::memory_buffer out;

  axom::fmt::format_to(std::back_inserter(out), "solid AxomMesh\n");
  for(auto t : tri_mesh.zoneSet)
  {
    auto tri = tri_mesh.znRel[t];
    const auto& v0 = tri_mesh.coords[tri[0]];
    const auto& v1 = tri_mesh.coords[tri[1]];
    const auto& v2 = tri_mesh.coords[tri[2]];
    const auto n = primal::Vector<double, 3>::cross_product(v1 - v0, v2 - v0);

    fmt::format_to(std::back_inserter(out), "  facet normal {} {} {}\n", n[0], n[1], n[2]);
    fmt::format_to(std::back_inserter(out), "    outer loop\n");
    fmt::format_to(std::back_inserter(out), "      vertex {} {} {}\n", v0[0], v0[1], v0[2]);
    fmt::format_to(std::back_inserter(out), "      vertex {} {} {}\n", v1[0], v1[1], v1[2]);
    fmt::format_to(std::back_inserter(out), "      vertex {} {} {}\n", v2[0], v2[1], v2[2]);
    fmt::format_to(std::back_inserter(out), "    endloop\n");
    fmt::format_to(std::back_inserter(out), "  endfacet\n");
  }
  fmt::format_to(std::back_inserter(out), "endsolid AxomMesh\n");

  std::ofstream stl_file(filename);
  if(!stl_file.is_open())
  {
    SLIC_ERROR("Failed to open STL file for writing: " << filename);
    return false;
  }
  stl_file.write(out.data(), out.size());
  SLIC_INFO("STL file written to: " << filename);

  return true;
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger(axom::slic::message::Debug);
  SLIC_INFO(axom::fmt::format("Axom Version: [{}]", axom::getVersion()));

  // Parse the command line arguments
  std::string input_file;
  std::string output_file;

  axom::CLI::App app {
    "mesh_converter -- converts from a volumetric ProE tetrahedral mesh to an STL mesh of "
    "the boundary triangles"};
  app.add_option("-i,--infile", input_file, "The input ProE file")
    ->required()
    ->check(axom::CLI::ExistingFile);
  app.add_option("-o,--outfile", output_file, "The output STL file")->required();

  double weld_threshold = 1e-6;
  app.add_option("-w,--weld-threshold", weld_threshold)
    ->description("Threshold for welding vertices")
    ->capture_default_str();

  bool verbose = false;
  app.add_flag("-v,--verbose", verbose, "Enable verbose output");

  CLI11_PARSE(app, argc, argv);

  // Load the ProE mesh, fixing tet orientations along the way
  TetMesh tet_mesh = loadProe(input_file, verbose);

  // Weld all vertices that are within \a weldThresole
  weldVertices(tet_mesh, weld_threshold, verbose);
  SLIC_INFO(axom::fmt::format("After welding, tet mesh has {} vertices and {} tets",
                              tet_mesh.nodeSet.size(),
                              tet_mesh.zoneSet.size()));

  // Create a tri mesh of (non-degenerate) boundary faces
  TriMesh tri_mesh;
  extractBoundaryFaces(tet_mesh, tri_mesh, verbose);
  SLIC_INFO(axom::fmt::format("Boundary triangle mesh has {} triangles", tri_mesh.zoneSet.size()));

  // Saves the result to an STL mesh
  bool success = writeSTL(tri_mesh, output_file);

  return success ? 0 : 1;
}