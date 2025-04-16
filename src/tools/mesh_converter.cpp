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

/**
 * This class represents a simplicial mesh
 * 
 * It uses slam to encode the vertices and simplices (zones) of a simplicial mesh 
 * as well as the connectivity relation from the simplices to vertices,
 * but does not represent other entities (e.g. edges, faces) 
 * or relations (like adjacency or vertex-simplex incidences).
 * 
 * It Can be used to represent both tetrahedral meshes (TOPO_DIM=3) 
 * and triangle meshes (TOPO_DIM=2).
 */
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

  // Define the SimplexType based on the topological dimension
  using SimplexType = std::conditional_t<
    TOPO_DIM == 3,
    primal::Tetrahedron<double, SPACE_DIM>,
    std::conditional_t<TOPO_DIM == 2, primal::Triangle<double, SPACE_DIM>, primal::Segment<double, SPACE_DIM>>>;

  static_assert(SPACE_DIM <= 3, "SPACE_DIM can be at most 3");
  static_assert(TOPO_DIM <= SPACE_DIM, "TOPO_DIM must be less than or equal to SPACE_DIM");
  static_assert(1 <= TOPO_DIM, "TOPO_DIM must be greater than 0");

public:
  SimplicialMesh(int num_verts = 0, int num_zones = 0) { reset(num_verts, num_zones); }

  /// Resets the internal sets and relations and ensures we have the right sizes
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

  /// Checks the validity of the mesh
  bool checkValid(bool verboseOutput = false) const
  {
    const bool nodes_valid = nodeSet.isValid(verboseOutput);
    const bool zones_valid = zoneSet.isValid(verboseOutput);
    const bool zn_valid =
      (nodeSet.size() > 0 && zoneSet.size() > 0) ? znRel.isValid(verboseOutput) : true;
    const bool coords_valid = coords.isValid(verboseOutput);

    return nodes_valid && zones_valid && zn_valid && coords_valid;
  }

  /**
   *  \brief Reindexes the mesh to remove duplicate vertices
   * 
   *  \param remap_verts a map from old vertex indices to new vertex indices
   *  \param retained_verts a bitset indicating which vertices to retain
   *  \param uniqueVertCount the number of unique vertices in the mesh
   */
  void reindexMesh(const axom::Array<PositionType>& remap_verts,
                   const slam::BitSet& retained_verts,
                   int uniqueVertCount)
  {
    // update node set and vertices
    {
      NodeSet newNodeSet(uniqueVertCount);
      Vertices newCoords(&newNodeSet);
      PositionType curr = 0;

      for(auto n = retained_verts.find_first(); n != slam::BitSet::npos;
          n = retained_verts.find_next(n))
      {
        newCoords[curr++] = coords[n];
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

  PointType getVertex(PositionType index) const { return coords[index]; }

  /// \brief Returns the segment corresponding to the given index (for TOPO_DIM == 1).
  template <int D = TOPO_DIM, typename std::enable_if<D == 1, int>::type = 0>
  SimplexType getSimplex(PositionType index) const
  {
    const auto zone = znRel[index];
    return SimplexType(coords[zone[0]], coords[zone[1]]);
  }

  /// \brief Returns the triangle corresponding to the given index (for TOPO_DIM == 2).
  template <int D = TOPO_DIM, typename std::enable_if<D == 2, int>::type = 0>
  SimplexType getSimplex(PositionType index) const
  {
    const auto zone = znRel[index];
    return SimplexType(coords[zone[0]], coords[zone[1]], coords[zone[2]]);
  }

  /// \brief Returns the tetrahedron corresponding to the given index (for TOPO_DIM == 3).
  template <int D = TOPO_DIM, typename std::enable_if<D == 3, int>::type = 0>
  SimplexType getSimplex(PositionType index) const
  {
    const auto zone = znRel[index];
    return SimplexType(coords[zone[0]], coords[zone[1]], coords[zone[2]], coords[zone[3]]);
  }

public:
  NodeSet nodeSet;
  ZoneSet zoneSet;

  ZoneToNodeRelation znRel;
  Vertices coords;

private:
  axom::Array<PositionType> m_zn;
};

using TetMesh = SimplicialMesh<3, 3>;  // Volumetric tetrahedral mesh in 3D
using TriMesh = SimplicialMesh<2, 3>;  // Surface triangle mesh in 3D

/// Loads the pro-e tet mesh and returns a TetMesh instance
/// Fixes orientations of the tets, but defers fixing degeneracies until the vertices have been welded
bool loadProe(const std::string& filename, TetMesh& slam_mesh, bool verbose)
{
  AXOM_ANNOTATE_SCOPE("load Pro-E file");

  namespace mint = axom::mint;

  SLIC_INFO(axom::fmt::format("Reading pro-e file: '{}' ...", filename));

  /// read pro-e mesh into a mint unstructured mesh
  auto mint_mesh = std::make_unique<mint::UnstructuredMesh<mint::SINGLE_SHAPE>>(3, axom::mint::TET);
  {
    axom::quest::ProEReader reader;
    reader.setFileName(filename);

    reader.read();
    SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                                "Input mesh has {:L} vertices and {:L} tets",
                                reader.getNumNodes(),
                                reader.getNumTets()));

    reader.getMesh(mint_mesh.get());
  }

  if(mint_mesh == nullptr)
  {
    SLIC_ERROR("Mesh is null after reading the proe file.");
    return false;
  }

  if(verbose)
  {
    std::string vtk_output_file = "proe_tet.vtk";
    mint::write_vtk(mint_mesh.get(), vtk_output_file);
    SLIC_INFO("Saved tet mesh to VTK file: " << vtk_output_file);
  }

  /// convert from mint to our slam-based TetMesh
  slam_mesh.reset(mint_mesh->getNumberOfNodes(), mint_mesh->getNumberOfCells());

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
    int orientation_fixes = 0;
    for(auto t : slam_mesh.zoneSet)
    {
      const axom::IndexType* connec = mint_mesh->getCellNodeIDs(t);
      auto tet = slam_mesh.znRel[t];
      tet[0] = connec[0];
      tet[1] = connec[1];
      tet[2] = connec[2];
      tet[3] = connec[3];

      // fix orientations -- tets should not have negative volumes
      if(slam_mesh.getSimplex(t).signedVolume() < 0.)
      {
        axom::utilities::swap(tet[2], tet[3]);
        ++orientation_fixes;
      }
    }
    SLIC_INFO_IF(verbose,
                 axom::fmt::format(axom::utilities::locale(),
                                   "Fixed orientations of {:L} tetrahedra",
                                   orientation_fixes));
  }

  slam_mesh.checkValid(verbose);

  return true;
}

template <typename MeshType>
void computeEdgeHistogram(const MeshType& mesh)
{
  AXOM_ANNOTATE_SCOPE(axom::fmt::format("compute {} edge histogram",
                                        std::is_same<MeshType, TriMesh>::value ? "tri" : "tet"));

  using MinMaxRange = primal::BoundingBox<double, 1>;
  using LengthType = MinMaxRange::PointType;

  MinMaxRange meshEdgeLenRange;

  // simple binning based on the exponent
  using LogHistogram = std::map<int, int>;
  LogHistogram edgeLenHist;  // Create histogram of edge lengths (log scale)

  using LogRangeMap = std::map<int, MinMaxRange>;
  LogRangeMap edgeLenRangeMap;  // Tracks range of edge lengths at each scale

  using Segment = primal::Segment<double, 3>;

  for(auto t : mesh.zoneSet)
  {
    auto simplex = mesh.znRel[t];

    // Add all edges of the simplex to the histogram
    // Note: this will count edges from neighboring simplexes more than once
    // but is still a useful way to undertand the distribution of edge lengths
    for(int i = 0; i < MeshType::NODES_PER_ZONE; ++i)
    {
      for(int j = i + 1; j < MeshType::NODES_PER_ZONE; ++j)
      {
        const double len = Segment(mesh.getVertex(simplex[i]), mesh.getVertex(simplex[j])).length();
        meshEdgeLenRange.addPoint(LengthType {len});

        int expBase2;
        std::frexp(len, &expBase2);
        edgeLenHist[expBase2]++;
        edgeLenRangeMap[expBase2].addPoint(LengthType {len});
      }
    }
  }

  // Print the edge length histogram
  {
    axom::fmt::memory_buffer out;
    axom::fmt::format_to(std::back_inserter(out), "Edge Length Histogram:\n");
    for(const auto& entry : edgeLenHist)
    {
      const auto key = entry.first;
      const auto count = entry.second;
      const auto& length_range = edgeLenRangeMap[key];
      axom::fmt::format_to(std::back_inserter(out),
                           axom::utilities::locale(),
                           "  2^{}: \t{:>10L} edges\t length range: [{:.4e} {:.4e}]\n",
                           key,
                           count,
                           length_range.getMin()[0],
                           length_range.getMax()[0]);
    }
    SLIC_INFO(axom::fmt::to_string(out));
  }
}

/// Welds all vertices in the tet mesh that are within \a weldThresh of each other
void weldVertices(TetMesh& mesh, double weld_thresh, bool verbose)
{
  AXOM_ANNOTATE_SCOPE("weld vertices");
  SLIC_ASSERT_MSG(weld_thresh > 0.,
                  "Welding threshold must be greater than 0. Passed in value was " << weld_thresh);

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
  for(const auto& offset : {TetMesh::PointType(0.), TetMesh::PointType(weld_thresh / 2.)})
  {
    // We're going to find the unique indices for the welded vertices
    const int num_verts = mesh.nodeSet.size();
    int unique_vert_count = 0;
    int weld_count = 0;

    // Set up the lattice for quantizing points to an integer lattice
    TetMesh::PointType origin(meshBB.getMin().array() - offset.array());
    Lattice3 lattice(origin, Lattice3::SpaceVector(TetMesh::PointType(weld_thresh)));

    // A hashmap from GridCells to the new vertex indices
    using GridCellToIndexMap = std::unordered_map<GridCell, IdxType, decltype(point_hash)>;
    GridCellToIndexMap vertexIndexMap(num_verts, point_hash);

    axom::Array<TetMesh::PositionType> vertex_remap(num_verts);
    slam::BitSet retained_verts(num_verts);
    for(auto n : mesh.nodeSet)
    {
      // find the new vertex index; if not already present, we'll keep this vertex
      auto res = vertexIndexMap.insert({lattice.gridCell(mesh.coords[n]), n});
      if(res.second == true)
      {
        vertex_remap[n] = unique_vert_count++;
        retained_verts.set(n);
      }
      else
      {
        const auto prev_vertex_idx = res.first->second;
        SLIC_ASSERT(prev_vertex_idx < n);
        weld_count++;
        vertex_remap[n] = vertex_remap[prev_vertex_idx];
      }
    }

    // Finally, update the mesh using the new vertex mapping
    mesh.reindexMesh(vertex_remap, retained_verts, unique_vert_count);

    SLIC_INFO_IF(verbose,
                 axom::fmt::format(axom::utilities::locale(),
                                   "Welded {:L} vertices while reindexing with offset = {}. Mesh "
                                   "now has {:L} vertices and {:L} tets.",
                                   weld_count,
                                   offset[0],
                                   mesh.nodeSet.size(),
                                   mesh.zoneSet.size()));
  }
}

/// Creates a \a tri_mesh out of the boundary faces of the \a tet_mesh
void extractBoundaryFaces(const TetMesh& tet_mesh, TriMesh& tri_mesh, bool verbose)
{
  AXOM_ANNOTATE_SCOPE("extract boundary faces");

  using PositionType = TetMesh::PositionType;
  using IndexTriple = std::tuple<PositionType, PositionType, PositionType>;

  const int max_facets = tet_mesh.zoneSet.size() * 4;
  axom::Array<std::pair<IndexTriple, int>> facets(0, max_facets);
  slam::BitSet unmatched(max_facets);

  auto is_degenerate = [](const auto& tet) -> bool {
    return tet[0] == tet[1] || tet[0] == tet[2] || tet[0] == tet[3]  //
      || tet[1] == tet[2] || tet[1] == tet[3]                        //
      || tet[2] == tet[3];
  };

  // we sort the vertex indices to avoid order dependence,
  // but keep the association to the original tet face so we don't mess up the orientations
  auto create_tuple = [](PositionType idx0, PositionType idx1, PositionType idx2) {
    int num_swaps = 0;
    if(idx0 > idx1)
    {
      axom::utilities::swap(idx0, idx1);
      ++num_swaps;
    }
    if(idx1 > idx2)
    {
      axom::utilities::swap(idx1, idx2);
      ++num_swaps;
    }
    if(idx0 > idx1)
    {
      axom::utilities::swap(idx0, idx1);
      ++num_swaps;
    }
    return std::make_pair(std::make_tuple(idx0, idx1, idx2), num_swaps);
  };

  {
    AXOM_ANNOTATE_SCOPE("find unmatched faces");

    int degenerate_count = 0;
    for(auto t : tet_mesh.zoneSet)
    {
      auto tet = tet_mesh.znRel[t];

      if(is_degenerate(tet))
      {
        ++degenerate_count;
        continue;
      }

      // assume input tet mesh is a manifold w/ boundaries
      // internal faces will appear twice; boundary faces once
      for(auto triple_face : {create_tuple(tet[1], tet[2], tet[3]),
                              create_tuple(tet[0], tet[3], tet[2]),
                              create_tuple(tet[0], tet[1], tet[3]),
                              create_tuple(tet[0], tet[2], tet[1])})
      {
        facets.emplace_back(triple_face);
      }
    }
    // Sort unmatched_faces using a custom sort function
    std::sort(facets.begin(),
              facets.end(),
              [](const std::pair<IndexTriple, int>& a, const std::pair<IndexTriple, int>& b) {
                return a.first < b.first;
              });

    auto differ = [](const auto& a, const auto& b) { return a.first != b.first; };

    for(int i = 1; i < facets.size() - 1; ++i)
    {
      if(differ(facets[i], facets[i - 1]) && differ(facets[i], facets[i + 1]))
      {
        unmatched.set(i);
      }
    }
    if(facets.size() > 2)
    {
      if(differ(facets[0], facets[1]))
      {
        unmatched.set(0);
      }
      if(differ(facets[facets.size() - 2], facets[facets.size() - 1]))
      {
        unmatched.set(facets.size() - 1);
      }
    }

    SLIC_INFO(axom::fmt::format(
      axom::utilities::locale(),
      "Retained {:L} boundary faces, out of {:L} (and skipped {:L} degenerate tets)",
      unmatched.count(),
      max_facets,
      degenerate_count));
  }

  /// Copy the boundary faces into the triangle mesh
  {
    AXOM_ANNOTATE_SCOPE("create triangle mesh");

    const int num_verts = tet_mesh.nodeSet.size();
    const int num_triangles = unmatched.count();
    tri_mesh.reset(num_verts, num_triangles);

    SLIC_INFO_IF(verbose,
                 axom::fmt::format(axom::utilities::locale(),
                                   "Triangle mesh has {:L} vertices and {:L} triangles",
                                   tri_mesh.nodeSet.size(),
                                   tri_mesh.zoneSet.size()));

    // 1. copy the vertices
    // Since STL doesn't care about vertex indices, our tri mesh will have all vertices from the tet mesh
    // Alternatively, we could process this a bit more and keep only the referenced boundary vertices
    for(auto n : tet_mesh.nodeSet)
    {
      tri_mesh.coords[n] = tet_mesh.coords[n];
    }

    // 2. copy the boundary triangles
    int curr = 0;
    for(auto idx = unmatched.find_first(); idx != slam::BitSet::npos; idx = unmatched.find_next(idx))
    {
      const auto& kv = facets[idx];
      const auto face_triple = kv.first;
      const bool swap_orient = kv.second % 2 == 0;

      // Note -- the triple was sorted, use the number of swaps to determine the orientation
      auto tri = tri_mesh.znRel[curr++];
      tri[0] = std::get<0>(face_triple);
      tri[swap_orient ? 1 : 2] = std::get<1>(face_triple);
      tri[swap_orient ? 2 : 1] = std::get<2>(face_triple);
    }
  }
}

/// Write the triangle mesh out as an ascii STL mesh
bool writeSTL(const TriMesh& tri_mesh, const std::string& filename)
{
  AXOM_ANNOTATE_SCOPE("write STL mesh");

  axom::fmt::memory_buffer out;

  axom::fmt::format_to(std::back_inserter(out), "solid AxomMesh\n");
  for(auto t : tri_mesh.zoneSet)
  {
    auto tri = tri_mesh.znRel[t];
    const auto& v0 = tri_mesh.coords[tri[0]];
    const auto& v1 = tri_mesh.coords[tri[1]];
    const auto& v2 = tri_mesh.coords[tri[2]];
    const auto n = primal::Vector<double, 3>::cross_product(v1 - v0, v2 - v0);

    axom::fmt::format_to(std::back_inserter(out), "  facet normal {} {} {}\n", n[0], n[1], n[2]);
    axom::fmt::format_to(std::back_inserter(out), "    outer loop\n");
    axom::fmt::format_to(std::back_inserter(out), "      vertex {} {} {}\n", v0[0], v0[1], v0[2]);
    axom::fmt::format_to(std::back_inserter(out), "      vertex {} {} {}\n", v1[0], v1[1], v1[2]);
    axom::fmt::format_to(std::back_inserter(out), "      vertex {} {} {}\n", v2[0], v2[1], v2[2]);
    axom::fmt::format_to(std::back_inserter(out), "    endloop\n");
    axom::fmt::format_to(std::back_inserter(out), "  endfacet\n");
  }
  axom::fmt::format_to(std::back_inserter(out), "endsolid AxomMesh\n");

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

  bool compute_edge_histogram = false;
  app.add_flag("--edge-histogram", compute_edge_histogram)
    ->description("Compute and display the edge histogram")
    ->capture_default_str();

  bool verbose = false;
  app.add_flag("-v,--verbose", verbose, "Enable verbose output");

  std::string annotationMode {"none"};
#ifdef AXOM_USE_CALIPER
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif

  CLI11_PARSE(app, argc, argv);

  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(annotationMode);

  AXOM_ANNOTATE_SCOPE("mesh_converter");

  // Load the ProE mesh, fixing tet orientations along the way
  TetMesh tet_mesh;
  bool load_success = loadProe(input_file, tet_mesh, verbose);

  if(compute_edge_histogram)
  {
    computeEdgeHistogram(tet_mesh);
  }

  // Weld all vertices that are within \a weldThresole
  weldVertices(tet_mesh, weld_threshold, verbose);
  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "After welding, tet mesh has {:L} vertices and {:L} tets",
                              tet_mesh.nodeSet.size(),
                              tet_mesh.zoneSet.size()));

  // Create a tri mesh of (non-degenerate) boundary faces
  TriMesh tri_mesh;
  extractBoundaryFaces(tet_mesh, tri_mesh, verbose);
  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Boundary triangle mesh has {:L} triangles",
                              tri_mesh.zoneSet.size()));

  if(compute_edge_histogram)
  {
    computeEdgeHistogram(tri_mesh);
  }

  // Saves the result to an STL mesh
  bool save_success = writeSTL(tri_mesh, output_file);

  return load_success && save_success ? 0 : 1;
}
