// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#ifdef AXOM_USE_MFEM

  #include "mfem.hpp"

  #include "gtest/gtest.h"

  #include "axom/sidre/core/sidre.hpp"
  #include "axom/sidre/core/MFEMSidreDataCollection.hpp"

using axom::sidre::Group;
using axom::sidre::MFEMSidreDataCollection;

const std::string COLL_NAME = "test_collection";
const double EPSILON = 1.0e-6;

TEST(sidre_datacollection, dc_alloc_no_mesh)
{
  MFEMSidreDataCollection sdc(COLL_NAME);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_owning_mesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc(COLL_NAME, &mesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_nonowning_mesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  bool owns_mesh = false;
  MFEMSidreDataCollection sdc(COLL_NAME, &mesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_register_empty_field)
{
  MFEMSidreDataCollection sdc(COLL_NAME);
  mfem::GridFunction gf;
  sdc.RegisterField("test_field", &gf);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_register_partial_field)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::H1_FECollection fec(2);
  mfem::FiniteElementSpace fes(&mesh, &fec);
  mfem::GridFunction gf(&fes);

  MFEMSidreDataCollection sdc(COLL_NAME, &mesh);
  sdc.RegisterField("test_field", &gf);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_update_state)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(COLL_NAME, &mesh);

  // Arbitrary values for the "state" part of blueprint
  sdc.SetCycle(3);
  sdc.SetTime(1.1);
  sdc.SetTimeStep(0.4);
  // Force refresh from the superclass
  sdc.UpdateStateToDS();

  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_save)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(COLL_NAME, &mesh);

  // The data produced isn't important for this, so just throw it in /tmp
  sdc.SetPrefixPath("/tmp/dc_save_test");
  sdc.Save();

  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_reload_gf)
{
  const std::string field_name = "test_field";
  // 2D mesh divided into triangles
  mfem::Mesh mesh(10, 10, mfem::Element::TRIANGLE);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(COLL_NAME, &mesh, owns_mesh);
  mfem::GridFunction gf_write(&fes, nullptr);

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, &gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write.ProjectCoefficient(three_and_a_half);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
  #endif

  sdc_writer.SetPrefixPath("/tmp/dc_reload_test");
  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  // No mesh is used here
  MFEMSidreDataCollection sdc_reader(COLL_NAME);

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
  #endif

  sdc_reader.SetPrefixPath("/tmp/dc_reload_test");
  sdc_reader.Load();

  // No need to reregister, it already exists
  auto gf_read = sdc_reader.GetField(field_name);

  // Make sure the gridfunction was actually read in
  EXPECT_LT(gf_read->ComputeL2Error(three_and_a_half), EPSILON);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_reload_mesh)
{
  const std::string field_name = "test_field";
  // 2D mesh divided into triangles
  mfem::Mesh mesh(10, 10, mfem::Element::TRIANGLE);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(COLL_NAME, &mesh, owns_mesh);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

  // Save some basic info about the mesh
  const int n_verts = sdc_writer.GetMesh()->GetNV();
  const int n_ele = sdc_writer.GetMesh()->GetNE();
  const int n_bdr_ele = sdc_writer.GetMesh()->GetNBE();

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
  #endif

  sdc_writer.SetPrefixPath("/tmp/dc_reload_test");
  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  // No mesh is used here
  MFEMSidreDataCollection sdc_reader(COLL_NAME);

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
  #endif

  sdc_reader.SetPrefixPath("/tmp/dc_reload_test");
  sdc_reader.Load();

  // Make sure the mesh was actually reconstructed
  EXPECT_EQ(sdc_reader.GetMesh()->GetNV(), n_verts);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNE(), n_ele);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNBE(), n_bdr_ele);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)

TEST(sidre_datacollection, dc_alloc_owning_parmesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc(COLL_NAME, &parmesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_nonowning_parmesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  bool owns_mesh = false;
  MFEMSidreDataCollection sdc(COLL_NAME, &parmesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

struct ParMeshGroupData
{
  int n_verts;
  int n_edges;
  int n_triangles;
  int n_quads;
  void check(const ParMeshGroupData& other)
  {
    EXPECT_EQ(n_verts, other.n_verts);
    EXPECT_EQ(n_edges, other.n_edges);
    EXPECT_EQ(n_triangles, other.n_triangles);
    EXPECT_EQ(n_quads, other.n_quads);
  }
};

static std::vector<ParMeshGroupData> getGroupData(const mfem::ParMesh& parmesh)
{
  std::vector<ParMeshGroupData> result(parmesh.GetNGroups());
  // remove when marked const in MFEM
  auto& non_const_parmesh = const_cast<mfem::ParMesh&>(parmesh);

  for(std::size_t i = 0; i < result.size(); i++)
  {
    result[i] = ParMeshGroupData {non_const_parmesh.GroupNVertices(i),
                                  non_const_parmesh.GroupNEdges(i),
                                  non_const_parmesh.GroupNTriangles(i),
                                  non_const_parmesh.GroupNQuadrilaterals(i)};
  }

  return result;
}

/**
 * @brief Helper method for testing that a parallel mesh is reconstructed correctly
 * @param [in] base_mesh The serial mesh object to distribute, save, and then reload
 * @param [in] part_method The partitioning method to use - in [0, 5]
 * @param [in] debug_print Whether to save pre- and post-reload mesh data to files for debugging
 */
static void testParallelMeshReload(mfem::Mesh& base_mesh,
                                   const int part_method = 1,
                                   bool debug_print = false)
{
  mfem::ParMesh parmesh(MPI_COMM_WORLD, base_mesh, nullptr, part_method);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int n_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);

  std::ofstream fout;
  if(debug_print)
  {
    fout.open("mesh_par_" + std::to_string(rank));
    parmesh.PrintSharedEntities("mesh_par");
    if(rank == 0)
    {
      // Prints local mesh on rank 0
      std::ofstream vtkstream("parmesh.vtk");
      parmesh.PrintVTK(vtkstream);
      std::ofstream part("partitioning.txt");
      auto partitioning = base_mesh.GeneratePartitioning(n_ranks, 5);
      for(int i = 0; i < base_mesh.GetNE(); i++)
      {
        part << "Element " << i << ": " << partitioning[i] << "\n";
      }
      delete[] partitioning;
    }
  }

  mfem::H1_FECollection fec(1, base_mesh.Dimension());
  mfem::ParFiniteElementSpace parfes(&parmesh, &fec);

  // The mesh must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  const bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(COLL_NAME, &parmesh, owns_mesh);
  if(debug_print && sdc_writer.GetBPGroup()->hasGroup("adjsets"))
  {
    sdc_writer.GetBPGroup()->print(fout);
  }

  // Save some basic info about the mesh
  const int n_verts = sdc_writer.GetMesh()->GetNV();
  const int n_ele = sdc_writer.GetMesh()->GetNE();
  const int n_bdr_ele = sdc_writer.GetMesh()->GetNBE();

  // ParMesh-specific info
  auto writer_pmesh = dynamic_cast<mfem::ParMesh*>(sdc_writer.GetMesh());
  ASSERT_NE(writer_pmesh, nullptr);
  const int n_groups = writer_pmesh->GetNGroups();
  const int n_face_neighbors = writer_pmesh->GetNFaceNeighbors();
  const int n_shared_faces = writer_pmesh->GetNSharedFaces();

  // Group-specific info
  auto writer_group_data = getGroupData(*writer_pmesh);

  sdc_writer.SetPrefixPath("/tmp/dc_par_reload_test");
  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  MFEMSidreDataCollection sdc_reader(COLL_NAME);

  // Needs to be set "manually" in order for everything to be loaded in properly
  sdc_reader.SetComm(MPI_COMM_WORLD);

  sdc_reader.SetPrefixPath("/tmp/dc_par_reload_test");
  sdc_reader.Load();

  // Make sure the mesh was actually reconstructed
  EXPECT_EQ(sdc_reader.GetMesh()->GetNV(), n_verts);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNE(), n_ele);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNBE(), n_bdr_ele);

  if(debug_print)
  {
    auto reader_pmesh = dynamic_cast<mfem::ParMesh*>(sdc_reader.GetMesh());
    if(reader_pmesh)
    {
      reader_pmesh->PrintSharedEntities("mesh_par_reload_");
    }
  }

  // Currently a mesh cannot be identified as a ParMesh in a vacuously
  // distributed setup
  if(n_ranks > 1)
  {
    // Make sure the ParMesh was reconstructed
    auto reader_pmesh = dynamic_cast<mfem::ParMesh*>(sdc_reader.GetMesh());
    ASSERT_NE(reader_pmesh, nullptr);
    EXPECT_EQ(reader_pmesh->GetNGroups(), n_groups);
    EXPECT_EQ(reader_pmesh->GetNFaceNeighbors(), n_face_neighbors);
    EXPECT_EQ(reader_pmesh->GetNSharedFaces(), n_shared_faces);

    auto reader_group_data = getGroupData(*reader_pmesh);
    EXPECT_EQ(writer_group_data.size(), reader_group_data.size());
    // std::equal would be nice here, but we would lose the gtest
    // diagnostics I think
    for(std::size_t i = 0; i < writer_group_data.size(); i++)
    {
      writer_group_data[i].check(reader_group_data[i]);
    }
  }

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

/**
 * @brief Wrapper for testing mesh reconstruction across all partitioning methods
 * @param [in] base_mesh The serial mesh object to distribute, save, and then reload
 * @param [in] debug_print Whether to save pre- and post-reload mesh data to files for debugging
 */
static void testParallelMeshReloadAllPartitionings(mfem::Mesh& base_mesh,
                                                   bool debug_print = false)
{
  for(int part_method = 0; part_method <= 5; part_method++)
  {
    testParallelMeshReload(base_mesh, part_method, debug_print);
  }
}

TEST(sidre_datacollection, dc_par_reload_gf)
{
  const std::string field_name = "test_field";
  // 3D tet mesh
  mfem::Mesh mesh(2, 2, 2, mfem::Element::TETRAHEDRON);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);

  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::ParFiniteElementSpace parfes(&parmesh, &fec);

  // The mesh must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(COLL_NAME, &parmesh, owns_mesh);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  mfem::ParGridFunction gf_write(&parfes, static_cast<double*>(nullptr));

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, &gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write.ProjectCoefficient(three_and_a_half);

  sdc_writer.SetPrefixPath("/tmp/dc_par_reload_test");
  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  MFEMSidreDataCollection sdc_reader(COLL_NAME);

  // Needs to be set "manually" in order for everything to be loaded in properly
  sdc_reader.SetComm(MPI_COMM_WORLD);

  sdc_reader.SetPrefixPath("/tmp/dc_par_reload_test");
  sdc_reader.Load();

  auto gf_read = sdc_reader.GetField(field_name);

  // Can only reconstruct if more than trivially parallel
  int n_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
  if(n_ranks > 1)
  {
    EXPECT_TRUE(dynamic_cast<mfem::ParGridFunction*>(gf_read));
  }

  // Make sure the gridfunction was actually read in
  EXPECT_LT(gf_read->ComputeL2Error(three_and_a_half), EPSILON);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

// MFEM's partitioning logic looks like it's buggy for 1D messages
// Probably not a problem since 1D meshes are rare

TEST(sidre_datacollection, dc_par_reload_mesh_2D_small)
{
  // 2D mesh divided into triangles
  mfem::Mesh mesh(10, 10, mfem::Element::TRIANGLE);
  testParallelMeshReloadAllPartitionings(mesh);
}

TEST(sidre_datacollection, dc_par_reload_mesh_2D_large)
{
  // 2D mesh divided into triangles
  mfem::Mesh mesh(100, 100, mfem::Element::TRIANGLE);
  testParallelMeshReloadAllPartitionings(mesh);
}

TEST(sidre_datacollection, dc_par_reload_mesh_3D_small_tet)
{
  // 3D mesh divided into tetrahedra
  mfem::Mesh mesh(2, 2, 2, mfem::Element::TETRAHEDRON);
  testParallelMeshReloadAllPartitionings(mesh);
}

TEST(sidre_datacollection, dc_par_reload_mesh_3D_medium_tet)
{
  // 3D mesh divided into tetrahedra
  mfem::Mesh mesh(10, 10, 10, mfem::Element::TETRAHEDRON);
  testParallelMeshReloadAllPartitionings(mesh);
}

TEST(sidre_datacollection, dc_par_reload_mesh_3D_small_hex)
{
  // 3D mesh divided into hexahedra
  mfem::Mesh mesh(2, 2, 2, mfem::Element::HEXAHEDRON);
  testParallelMeshReloadAllPartitionings(mesh);
}

TEST(sidre_datacollection, dc_par_reload_mesh_3D_medium_hex)
{
  // 3D mesh divided into hexahedra
  mfem::Mesh mesh(10, 10, 10, mfem::Element::HEXAHEDRON);
  testParallelMeshReloadAllPartitionings(mesh);
}

    //----------------------------------------------------------------------
    #include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}

  #endif

#endif
