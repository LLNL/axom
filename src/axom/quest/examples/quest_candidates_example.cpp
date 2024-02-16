// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

//-----------------------------------------------------------------------------
///
/// file: quest_candidates_examples.cpp
///
/// This example takes as input two Blueprint unstructured hex meshes, and
/// finds the candidates of intersection between the meshes using a
/// spatial index, either a Bounding Volume Hierarchy or an Implicit Grid.
/// The example supports HIP and CUDA execution through RAJA.
//-----------------------------------------------------------------------------

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"

#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"
#else
  #error This example requires axom to be configured with RAJA support
#endif

#ifdef AXOM_USE_UMPIRE
  #include "umpire/Umpire.hpp"
#else
  #error This example requires axom to be configured with Umpire support
#endif

#ifdef AXOM_USE_CONDUIT
  #include "conduit_relay.hpp"
  #include "conduit_blueprint.hpp"
#else
  #error This example requires axom to be configured with Conduit support
#endif

#include "axom/mint.hpp"
#include "axom/primal.hpp"
#include "axom/spin.hpp"
#include "axom/quest.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <memory>
#include <array>

using seq_exec = axom::SEQ_EXEC;

using UMesh = axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>;
using IndexPair = std::pair<axom::IndexType, axom::IndexType>;

// clang-format off
#if defined(AXOM_USE_OPENMP)
  using omp_exec = axom::OMP_EXEC;
#else
  using omp_exec = seq_exec;
#endif

#if defined(AXOM_USE_CUDA)
  constexpr int CUDA_BLK_SZ = 256;
  using cuda_exec = axom::CUDA_EXEC<CUDA_BLK_SZ>;
#else
  using cuda_exec = seq_exec;
#endif

// Hip exec included from IntersectionShaper header

// clang-format on

//-----------------------------------------------------------------------------
/// Basic RAII utility class for initializing and finalizing slic logger
//-----------------------------------------------------------------------------
struct BasicLogger
{
  BasicLogger()
  {
    namespace slic = axom::slic;

    // Initialize the SLIC logger
    slic::initialize();
    slic::setLoggingMsgLevel(slic::message::Debug);

    // Customize logging levels and formatting
    const std::string slicFormatStr = "[<LEVEL>] <MESSAGE> \n";

    slic::addStreamToMsgLevel(new slic::GenericOutputStream(&std::cerr),
                              slic::message::Error);
    slic::addStreamToMsgLevel(
      new slic::GenericOutputStream(&std::cerr, slicFormatStr),
      slic::message::Warning);

    auto* compactStream =
      new slic::GenericOutputStream(&std::cout, slicFormatStr);
    slic::addStreamToMsgLevel(compactStream, slic::message::Info);
    slic::addStreamToMsgLevel(compactStream, slic::message::Debug);
  }

  ~BasicLogger() { axom::slic::finalize(); }
};

//-----------------------------------------------------------------------------
/// Struct to help with parsing and storing command line args
//-----------------------------------------------------------------------------
enum class RuntimePolicy
{
  raja_seq = 1,
  raja_omp = 2,
  raja_cuda = 3,
  raja_hip = 4
};

struct Input
{
  static const std::set<std::string> s_validMethods;
  static const std::map<std::string, RuntimePolicy> s_validPolicies;

  std::string mesh_file_first {""};
  std::string mesh_file_second {""};
  std::string method {"bvh"};
  int resolution {0};
  bool verboseOutput {false};
  RuntimePolicy policy {RuntimePolicy::raja_seq};

  void parse(int argc, char** argv, axom::CLI::App& app);
  bool isVerbose() const { return verboseOutput; }
};

void Input::parse(int argc, char** argv, axom::CLI::App& app)
{
  app.add_option("-i, --infile", mesh_file_first)
    ->description(
      "The first input Blueprint mesh file to insert into spatial index")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app.add_option("-q, --queryfile", mesh_file_second)
    ->description("The second input Blueprint mesh file to query spatial index")
    ->required()
    ->check(axom::CLI::ExistingFile);

  app
    .add_option("-r,--resolution",
                resolution,
                "With '-m implicit', set resolution of implicit grid. \n"
                "Set to less than 1 to use the implicit grid spatial index\n"
                "with a resolution of the cube root of the number of\n"
                "hexes.")
    ->capture_default_str();

  app.add_flag("-v,--verbose", verboseOutput)
    ->description("Increase logging verbosity")
    ->capture_default_str();

  app.add_option("-p, --policy", policy)
    ->description(
      "Execution policy."
      "\nSet to 'raja_seq' or 1 to use the RAJA sequential policy."
#ifdef AXOM_USE_OPENMP
      "\nSet to 'raja_omp' or 2 to use the RAJA openmp policy."
#endif
#ifdef AXOM_USE_CUDA
      "\nSet to 'raja_cuda' or 3 to use the RAJA cuda policy."
#endif
#ifdef AXOM_USE_HIP
      "\nSet to 'raja_hip' or 4 to use the RAJA hip policy."
#endif
      )
    ->capture_default_str()
    ->transform(axom::CLI::CheckedTransformer(Input::s_validPolicies));

  app
    .add_option(
      "-m, --method",
      method,
      "Method to use. \n"
      "Set to 'bvh' to use the bounding volume hierarchy spatial index.\n"
      "Set to 'implicit' to use the implicit grid spatial index.")
    ->capture_default_str()
    ->check(axom::CLI::IsMember {Input::s_validMethods});

  app.get_formatter()->column_width(40);

  app.parse(argc, argv);  // Could throw an exception

  // Output parsed information
  SLIC_INFO(axom::fmt::format(
    R"(
     Parsed parameters:
      * First Blueprint mesh to insert: '{}'
      * Second Blueprint mesh to query: '{}'
      * Verbose logging: {}
      * Spatial method: '{}'
      * Resolution: '{}'
      * Runtime execution policy: '{}'
      )",
    mesh_file_first,
    mesh_file_second,
    verboseOutput,
    method == "bvh" ? "Bounding Volume Hierarchy (BVH)" : "Implicit Grid",
    method == "bvh" ? "Not Applicable" : std::to_string(resolution),
    policy == RuntimePolicy::raja_omp ? "raja_omp"
                                      : (policy == RuntimePolicy::raja_cuda)
        ? "raja_cuda"
        : (policy == RuntimePolicy::raja_hip) ? "raja_hip" : "raja_seq"));
}

const std::map<std::string, RuntimePolicy> Input::s_validPolicies(
  {{"raja_seq", RuntimePolicy::raja_seq}
#ifdef AXOM_USE_OPENMP
   ,
   {"raja_omp", RuntimePolicy::raja_omp}
#endif
#ifdef AXOM_USE_CUDA
   ,
   {"raja_cuda", RuntimePolicy::raja_cuda}
#endif
#ifdef AXOM_USE_HIP
   ,
   {"raja_hip", RuntimePolicy::raja_hip}
#endif
  });

const std::set<std::string> Input::s_validMethods({
  "bvh",
  "implicit",
});

//-----------------------------------------------------------------------------
/// Basic hexahedron mesh to be used in our application
//-----------------------------------------------------------------------------
struct HexMesh
{
  using Point = axom::primal::Point<double, 3>;
  using Hexahedron = axom::primal::Hexahedron<double, 3>;
  using BoundingBox = axom::primal::BoundingBox<double, 3>;

  axom::IndexType numHexes() const { return m_hexes.size(); }
  axom::Array<Hexahedron>& hexes() { return m_hexes; }
  const axom::Array<Hexahedron>& hexes() const { return m_hexes; }

  BoundingBox& meshBoundingBox() { return m_meshBoundingBox; }
  const BoundingBox& meshBoundingBox() const { return m_meshBoundingBox; }

  axom::Array<BoundingBox>& hexBoundingBoxes() { return m_hexBoundingBoxes; }
  const axom::Array<BoundingBox>& hexBoundingBoxes() const
  {
    return m_hexBoundingBoxes;
  }

  axom::Array<Hexahedron> m_hexes;
  axom::Array<BoundingBox> m_hexBoundingBoxes;
  BoundingBox m_meshBoundingBox;
};

HexMesh loadBlueprintHexMesh(const std::string& mesh_path,
                             bool verboseOutput = false)
{
  HexMesh hexMesh;

  // Load Blueprint mesh into Conduit node
  conduit::Node n_load;
  conduit::relay::io::blueprint::read_mesh(mesh_path, n_load);

  // Verify this is a hexahedral mesh
  std::string shape = n_load[0]["topologies/topo/elements/shape"].as_string();
  if(shape != "hex")
  {
    SLIC_ERROR("A hex mesh was expected!");
  }

  const int HEX_OFFSET = 8;

  int num_nodes =
    (n_load[0]["coordsets/coords/values/x"]).dtype().number_of_elements();

  int connectivity_size =
    (n_load[0]["topologies/topo/elements/connectivity"]).dtype().number_of_elements();

  // Sanity check for number of cells
  int cell_calc_from_nodes =
    std::round(std::pow(std::pow(num_nodes, 1.0 / 3.0) - 1, 3));
  int cell_calc_from_connectivity = connectivity_size / HEX_OFFSET;
  if(cell_calc_from_nodes != cell_calc_from_connectivity)
  {
    SLIC_ERROR("Number of cells is not expected!\n"
               << "First calculation is " << cell_calc_from_nodes
               << " and second calculation is " << cell_calc_from_connectivity);
  }

  // extract hexes into an axom::Array
  int* connectivity = n_load[0]["topologies/topo/elements/connectivity"].value();

  double* x_vals = n_load[0]["coordsets/coords/values/x"].value();
  double* y_vals = n_load[0]["coordsets/coords/values/y"].value();
  double* z_vals = n_load[0]["coordsets/coords/values/z"].value();

  const int numCells = connectivity_size / HEX_OFFSET;
  hexMesh.m_hexes.reserve(numCells);
  HexMesh::Hexahedron hex;
  axom::Array<HexMesh::Point> hexPoints(HEX_OFFSET);

  for(int i = 0; i < numCells; ++i)
  {
    for(int j = 0; j < HEX_OFFSET; j++)
    {
      int offset = i * HEX_OFFSET;
      hexPoints[j] = HexMesh::Point({x_vals[connectivity[offset + j]],
                                     y_vals[connectivity[offset + j]],
                                     z_vals[connectivity[offset + j]]});
    }
    hex = HexMesh::Hexahedron(hexPoints);
    hexMesh.m_hexes.emplace_back(hex);
  }

  // compute and store hex bounding boxes and mesh bounding box
  hexMesh.m_hexBoundingBoxes.reserve(numCells);
  for(const auto& hex : hexMesh.hexes())
  {
    hexMesh.m_hexBoundingBoxes.emplace_back(
      axom::primal::compute_bounding_box(hex));
    hexMesh.m_meshBoundingBox.addBox(hexMesh.m_hexBoundingBoxes.back());
  }

  SLIC_INFO(
    axom::fmt::format("Mesh bounding box is {}.\n", hexMesh.meshBoundingBox()));

  // Optional verbose output that writes Blueprint mesh to vtk
  if(verboseOutput)
  {
    UMesh* mesh = new UMesh(3, axom::mint::HEX);

    // Append mesh nodes
    for(int i = 0; i < num_nodes; i++)
    {
      mesh->appendNode(x_vals[i], y_vals[i], z_vals[i]);
    }

    // Append mesh cells
    for(int i = 0; i < numCells; i++)
    {
      const axom::IndexType cell[] = {
        connectivity[i * HEX_OFFSET],
        connectivity[(i * HEX_OFFSET) + 1],
        connectivity[(i * HEX_OFFSET) + 2],
        connectivity[(i * HEX_OFFSET) + 3],
        connectivity[(i * HEX_OFFSET) + 4],
        connectivity[(i * HEX_OFFSET) + 5],
        connectivity[(i * HEX_OFFSET) + 6],
        connectivity[(i * HEX_OFFSET) + 7],
      };

      mesh->appendCell(cell);
    }

    // Write out to vtk for test viewing
    SLIC_INFO("Writing out Blueprint mesh to test.vtk for debugging...");
    axom::utilities::Timer timer(true);
    axom::mint::write_vtk(mesh, "test.vtk");
    timer.stop();
    SLIC_INFO(axom::fmt::format(
      "Writing out Blueprint mesh to test.vtk took {:4.3} seconds.",
      timer.elapsedTimeInSec()));

    delete mesh;
    mesh = nullptr;
  }  // end of verbose output

  return hexMesh;
}  // end of loadBlueprintHexMesh

template <typename ExecSpace>
axom::Array<IndexPair> findCandidatesBVH(const HexMesh& insertMesh,
                                         const HexMesh& queryMesh)
{
  SLIC_INFO("Running BVH candidates algorithm in execution Space: "
            << axom::execution_space<ExecSpace>::name());

  using HexArray = axom::Array<typename HexMesh::Hexahedron>;
  using BBoxArray = axom::Array<typename HexMesh::BoundingBox>;
  using IndexArray = axom::Array<axom::IndexType>;
  constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();

  axom::Array<IndexPair> candidatePairs;

  // Get ids of necessary allocators
  const int host_allocator =
    axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  const int kernel_allocator = on_device
    ? axom::getUmpireResourceAllocatorID(umpire::resource::Device)
    : axom::execution_space<ExecSpace>::allocatorID();

  // Copy the insert-BVH hexes to the device, if necessary
  // Either way, insert_hexes_v will be a view w/ data in the correct space
  auto& insert_hexes_h = insertMesh.hexes();
  HexArray insert_hexes_d =
    on_device ? HexArray(insert_hexes_h, kernel_allocator) : HexArray();

  // Copy the insert-BVH bboxes to the device, if necessary
  // Either way, insert_bbox_v will be a view w/ data in the correct space
  auto& insert_bbox_h = insertMesh.hexBoundingBoxes();
  BBoxArray insert_bbox_d =
    on_device ? BBoxArray(insert_bbox_h, kernel_allocator) : BBoxArray();
  auto insert_bbox_v = on_device ? insert_bbox_d.view() : insert_bbox_h.view();

  // Copy the query-BVH hexes to the device, if necessary
  // Either way, query_hexes_v will be a view w/ data in the correct space
  auto& query_hexes_h = queryMesh.hexes();
  HexArray query_hexes_d =
    on_device ? HexArray(query_hexes_h, kernel_allocator) : HexArray();

  // Copy the query-BVH bboxes to the device, if necessary
  // Either way, bbox_v will be a view w/ data in the correct space
  auto& query_bbox_h = queryMesh.hexBoundingBoxes();
  BBoxArray query_bbox_d =
    on_device ? BBoxArray(query_bbox_h, kernel_allocator) : BBoxArray();
  auto query_bbox_v = on_device ? query_bbox_d.view() : query_bbox_h.view();

  axom::utilities::Timer timer;

  // Initialize a BVH tree over the insert mesh bounding boxes
  timer.start();
  axom::spin::BVH<3, ExecSpace, double> bvh;
  bvh.setAllocatorID(kernel_allocator);
  bvh.initialize(insert_bbox_v, insert_bbox_v.size());
  timer.stop();
  SLIC_INFO(axom::fmt::format("0: Initializing BVH took {:4.3} seconds.",
                              timer.elapsedTimeInSec()));

  // Search for candidate bounding boxes of hexes to query;
  timer.start();
  IndexArray offsets_d(query_bbox_v.size(), query_bbox_v.size(), kernel_allocator);
  IndexArray counts_d(query_bbox_v.size(), query_bbox_v.size(), kernel_allocator);
  IndexArray candidates_d(0, 0, kernel_allocator);

  auto offsets_v = offsets_d.view();
  auto counts_v = counts_d.view();
  bvh.findBoundingBoxes(offsets_v,
                        counts_v,
                        candidates_d,
                        query_bbox_v.size(),
                        query_bbox_v);

  timer.stop();
  SLIC_INFO(axom::fmt::format(
    "1: Querying candidate bounding boxes took {:4.3} seconds.",
    timer.elapsedTimeInSec()));

  // Initialize candidatePairs to return
  timer.start();

  IndexArray offsets_h(offsets_d, host_allocator);
  IndexArray counts_h(counts_d, host_allocator);
  IndexArray candidates_h(candidates_d, host_allocator);

  for(int i = 0; i < queryMesh.numHexes(); i++)
  {
    for(int j = 0; j < counts_h[i]; j++)
    {
      candidatePairs.emplace_back(
        std::make_pair(i, candidates_h[offsets_h[i] + j]));
    }
  }
  timer.stop();

  SLIC_INFO(axom::fmt::format(
    "2: Initializing candidate pairs on host took {:4.3} seconds.",
    timer.elapsedTimeInSec()));

  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              R"(Stats for query
    -- Number of insert-BVH mesh hexes {:L}
    -- Number of query mesh hexes {:L}
    -- Total possible candidates {:L}
    -- Candidates from BVH query {:L}
    )",
                              insertMesh.numHexes(),
                              queryMesh.numHexes(),
                              1.0 * insertMesh.numHexes() * queryMesh.numHexes(),
                              candidates_h.size()));

  return candidatePairs;
}  // end of findCandidatesBVH for Blueprint Meshes

template <typename ExecSpace>
axom::Array<IndexPair> findCandidatesImplicit(const HexMesh& insertMesh,
                                              const HexMesh& queryMesh,
                                              int resolution)
{
  axom::Array<IndexPair> candidatePairs;

  SLIC_INFO("Running Implicit Grid candidates algorithm in execution Space: "
            << axom::execution_space<ExecSpace>::name());

  using HexArray = axom::Array<typename HexMesh::Hexahedron>;
  using BBoxArray = axom::Array<typename HexMesh::BoundingBox>;
  using IndexArray = axom::Array<int>;
  constexpr bool on_device = axom::execution_space<ExecSpace>::onDevice();

  // Get ids of necessary allocators
  const int host_allocator =
    axom::getUmpireResourceAllocatorID(umpire::resource::Host);
  const int kernel_allocator = on_device
    ? axom::getUmpireResourceAllocatorID(umpire::resource::Device)
    : axom::execution_space<ExecSpace>::allocatorID();

  // Copy the insert hexes to the device, if necessary
  // Either way, insert_hexes_v will be a view w/ data in the correct space
  auto& insert_hexes_h = insertMesh.hexes();
  HexArray insert_hexes_d =
    on_device ? HexArray(insert_hexes_h, kernel_allocator) : HexArray();

  // Copy the insert bboxes to the device, if necessary
  // Either way, insert_bbox_v will be a view w/ data in the correct space
  auto& insert_bbox_h = insertMesh.hexBoundingBoxes();
  BBoxArray insert_bbox_d =
    on_device ? BBoxArray(insert_bbox_h, kernel_allocator) : BBoxArray();
  auto insert_bbox_v = on_device ? insert_bbox_d.view() : insert_bbox_h.view();

  // Bounding box of entire insert mesh
  HexMesh::BoundingBox insert_mesh_bbox_h = insertMesh.meshBoundingBox();

  // Copy the query hexes to the device, if necessary
  // Either way, query_hexes_v will be a view w/ data in the correct space
  auto& query_hexes_h = queryMesh.hexes();
  HexArray query_hexes_d =
    on_device ? HexArray(query_hexes_h, kernel_allocator) : HexArray();

  // Copy the query bboxes to the device, if necessary
  // Either way, bbox_v will be a view w/ data in the correct space
  auto& query_bbox_h = queryMesh.hexBoundingBoxes();
  BBoxArray query_bbox_d =
    on_device ? BBoxArray(query_bbox_h, kernel_allocator) : BBoxArray();
  auto query_bbox_v = on_device ? query_bbox_d.view() : query_bbox_h.view();

  axom::utilities::Timer timer;
  timer.start();

  // If given resolution is less than one, use the cube root of the
  // number of hexes
  if(resolution < 1)
  {
    resolution = (int)(1 + std::pow(insertMesh.numHexes(), 1 / 3.));
  }

  const axom::primal::Point<int, 3> resolutions(resolution);

  axom::spin::ImplicitGrid<3, ExecSpace, int> gridIndex(insert_mesh_bbox_h,
                                                        &resolutions,
                                                        insertMesh.numHexes(),
                                                        kernel_allocator);
  gridIndex.insert(insertMesh.numHexes(), insert_bbox_v.data());
  timer.stop();
  SLIC_INFO(
    axom::fmt::format("0: Initializing Implicit Grid took {:4.3} seconds.",
                      timer.elapsedTimeInSec()));

  timer.start();
  IndexArray offsets_d(query_bbox_v.size(), query_bbox_v.size(), kernel_allocator);
  IndexArray counts_d(query_bbox_v.size(), query_bbox_v.size(), kernel_allocator);
  IndexArray candidates_d(0, 0, kernel_allocator);

  auto offsets_v = offsets_d.view();
  auto counts_v = counts_d.view();

  // First pass: get number of bounding box candidates for each query bounding box
  // Logic here mirrors the BVH two-pass query
  const auto grid_device = gridIndex.getQueryObject();
  using reduce_pol = typename axom::execution_space<ExecSpace>::reduce_policy;
  RAJA::ReduceSum<reduce_pol, int> totalCandidatePairs(0);

  axom::for_all<ExecSpace>(
    queryMesh.numHexes(),
    AXOM_LAMBDA(int icell) {
      int count = 0;

      // Define a function that is called on every candidate reached during
      // traversal. The function below simply counts the number of candidates
      // that intersect with the given query bounding box.
      auto isBBIntersect = [&](int candidateIdx) {
        if(axom::primal::intersect(insert_bbox_v[candidateIdx],
                                   query_bbox_v[icell]))
        {
          count++;
        }
      };

      // Call visitCandidates to iterate through the candidates
      grid_device.visitCandidates(query_bbox_v[icell], isBBIntersect);

      // Store the number of intersections for each query bounding box
      counts_v[icell] = count;
      totalCandidatePairs += count;
    });

  // Generate offsets from the counts
  using exec_pol = typename axom::execution_space<ExecSpace>::loop_policy;
  RAJA::exclusive_scan<exec_pol>(
    RAJA::make_span(counts_v.data(), queryMesh.numHexes()),
    RAJA::make_span(offsets_v.data(), queryMesh.numHexes()),
    RAJA::operators::plus<int> {});

  timer.stop();
  SLIC_INFO(axom::fmt::format(
    "1: Counting candidate bounding boxes took {:4.3} seconds.",
    timer.elapsedTimeInSec()));

  // Initialize candidatePairs to return
  timer.start();

  // Allocate arrays for intersection pairs
  IndexArray firstPair_d(totalCandidatePairs.get(),
                         totalCandidatePairs.get(),
                         kernel_allocator);
  IndexArray secondPair_d(totalCandidatePairs.get(),
                          totalCandidatePairs.get(),
                          kernel_allocator);
  auto first_pair_v = firstPair_d.view();
  auto second_pair_v = secondPair_d.view();

  // Second pass: fill candidates array pairs on device
  axom::for_all<ExecSpace>(
    queryMesh.numHexes(),
    AXOM_LAMBDA(axom::IndexType icell) {
      axom::IndexType offset = offsets_v[icell];

      // Store the intersection candidate
      auto fillCandidates = [&](int candidateIdx) {
        if(axom::primal::intersect(insert_bbox_v[candidateIdx],
                                   query_bbox_v[icell]))
        {
          first_pair_v[offset] = icell;
          second_pair_v[offset] = candidateIdx;
          offset++;
        }
      };

      grid_device.visitCandidates(query_bbox_v[icell], fillCandidates);
    });

  timer.stop();

  SLIC_INFO(axom::fmt::format(
    "2: Initializing candidate pairs on device took {:4.3} seconds.",
    timer.elapsedTimeInSec()));

  // copy results back to host and into return vector
  timer.start();

  IndexArray candidates_h[2] = {
    on_device ? IndexArray(firstPair_d, host_allocator) : IndexArray(),
    on_device ? IndexArray(secondPair_d, host_allocator) : IndexArray()};

  auto candidate1_h_v = on_device ? candidates_h[0].view() : firstPair_d.view();
  auto candidate2_h_v = on_device ? candidates_h[1].view() : secondPair_d.view();

  for(int idx = 0; idx < totalCandidatePairs; ++idx)
  {
    candidatePairs.emplace_back(
      std::make_pair(candidate1_h_v[idx], candidate2_h_v[idx]));
  }

  timer.stop();

  SLIC_INFO(
    axom::fmt::format("3: Moving candidate pairs to host took {:4.3} seconds.",
                      timer.elapsedTimeInSec()));

  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              R"(Stats for query
    -- Number of insert mesh hexes {:L}
    -- Number of query mesh hexes {:L}
    -- Total possible candidates {:L}
    -- Candidates from Implicit Grid query {:L}
    )",
                              insertMesh.numHexes(),
                              queryMesh.numHexes(),
                              1.0 * insertMesh.numHexes() * queryMesh.numHexes(),
                              candidatePairs.size()));

  return candidatePairs;
}

int main(int argc, char** argv)
{
  // Initialize logger; use RAII so it will finalize at the end of the application
  BasicLogger logger;

  // Parse the command line arguments
  Input params;
  {
    axom::CLI::App app {"Blueprint Hex BVH mesh candidate tester"};
    try
    {
      params.parse(argc, argv, app);
    }
    catch(const axom::CLI::ParseError& e)
    {
      return app.exit(e);
    }
  }

  axom::utilities::Timer timer(true);

  // Update the logging level based on verbosity flag
  axom::slic::setLoggingMsgLevel(params.isVerbose() ? axom::slic::message::Debug
                                                    : axom::slic::message::Info);

  // Load Blueprint mesh to insert into spatial index
  SLIC_INFO(axom::fmt::format("Reading Blueprint file to insert: '{}'...\n",
                              params.mesh_file_first));

  HexMesh insert_mesh =
    loadBlueprintHexMesh(params.mesh_file_first, params.isVerbose());

  // Load Blueprint mesh for querying spatial index
  SLIC_INFO(axom::fmt::format("Reading Blueprint file to query: '{}'...\n",
                              params.mesh_file_second));

  HexMesh query_mesh =
    loadBlueprintHexMesh(params.mesh_file_second, params.isVerbose());

  timer.stop();

  SLIC_INFO(axom::fmt::format("Reading in Blueprint files took {:4.3} seconds.",
                              timer.elapsedTimeInSec()));

  // Check for candidates; results are returned as an array of index pairs
  axom::Array<IndexPair> candidatePairs;
  timer.start();

  if(params.method == "bvh")
  {
    switch(params.policy)
    {
    case RuntimePolicy::raja_omp:
#ifdef AXOM_USE_OPENMP
      candidatePairs = findCandidatesBVH<omp_exec>(insert_mesh, query_mesh);
#endif
      break;
    case RuntimePolicy::raja_cuda:
#ifdef AXOM_USE_CUDA
      candidatePairs = findCandidatesBVH<cuda_exec>(insert_mesh, query_mesh);
#endif
      break;
    case RuntimePolicy::raja_hip:
#ifdef AXOM_USE_HIP
      candidatePairs = findCandidatesBVH<hip_exec>(insert_mesh, query_mesh);
#endif
      break;
    default:  // RuntimePolicy::raja_seq
      candidatePairs = findCandidatesBVH<seq_exec>(insert_mesh, query_mesh);
      break;
    }
  }
  // Implicit Grid
  else
  {
    switch(params.policy)
    {
    case RuntimePolicy::raja_omp:
#ifdef AXOM_USE_OPENMP
      candidatePairs = findCandidatesImplicit<omp_exec>(insert_mesh,
                                                        query_mesh,
                                                        params.resolution);
#endif
      break;
    case RuntimePolicy::raja_cuda:
#ifdef AXOM_USE_CUDA
      candidatePairs = findCandidatesImplicit<cuda_exec>(insert_mesh,
                                                         query_mesh,
                                                         params.resolution);
#endif
      break;
    case RuntimePolicy::raja_hip:
#ifdef AXOM_USE_HIP
      candidatePairs = findCandidatesImplicit<hip_exec>(insert_mesh,
                                                        query_mesh,
                                                        params.resolution);
#endif
      break;
    default:  // RuntimePolicy::raja_seq
      candidatePairs = findCandidatesImplicit<seq_exec>(insert_mesh,
                                                        query_mesh,
                                                        params.resolution);
      break;
    }
  }
  timer.stop();

  SLIC_INFO(axom::fmt::format("Computing candidates took {:4.3} seconds.",
                              timer.elapsedTimeInSec()));
  SLIC_INFO(axom::fmt::format(axom::utilities::locale(),
                              "Mesh had {:L} candidates pairs",
                              candidatePairs.size()));

  // print first few pairs
  const int numCandidates = candidatePairs.size();
  if(numCandidates > 0 && params.isVerbose())
  {
    constexpr int MAX_PRINT = 20;
    if(numCandidates > MAX_PRINT)
    {
      candidatePairs.resize(MAX_PRINT);
      SLIC_INFO(axom::fmt::format("First {} candidate pairs: {} ...\n",
                                  MAX_PRINT,
                                  axom::fmt::join(candidatePairs, ", ")));
    }
    else
    {
      SLIC_INFO(axom::fmt::format("Candidate pairs: {}\n",
                                  axom::fmt::join(candidatePairs, ", ")));
    }

    // Write out candidate pairs
    SLIC_INFO("Writing out candidate pairs...");
    std::ofstream outf("candidates.txt");

    outf << candidatePairs.size() << " candidate pairs:" << std::endl;
    for(int i = 0; i < candidatePairs.size(); ++i)
    {
      outf << candidatePairs[i].first << " " << candidatePairs[i].second
           << std::endl;
    }
  }

  return 0;
}
