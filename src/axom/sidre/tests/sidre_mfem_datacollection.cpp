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
/*
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
*/

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
/*
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
*/

TEST(sidre_datacollection, dc_par_reload_gf)
{
  const std::string field_name = "test_field";
  // 1D mesh divided into 10 segments
  // mfem::Mesh mesh(10);
  // 2D mesh divided into triangles
  // mfem::Mesh mesh(10, 10, mfem::Element::TRIANGLE);
  mfem::Mesh mesh(2, 2, 2, mfem::Element::TETRAHEDRON);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::ofstream fout("mesh_par_" + std::to_string(rank));
  parmesh.ParPrint(fout);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::ParFiniteElementSpace parfes(&parmesh, &fec);

  // The mesh must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(COLL_NAME, &parmesh, owns_mesh);
  sdc_writer.GetBPGroup()->getGroup("adjsets/mesh/groups")->print(fout);

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

  std::ofstream fout_rel("mesh_par_reload_" + std::to_string(rank));
  dynamic_cast<mfem::ParMesh*>(sdc_reader.GetMesh())->ParPrint(fout_rel);

  auto gf_read = sdc_reader.GetField(field_name);

  // Make sure the gridfunction was actually read in
  EXPECT_LT(gf_read->ComputeL2Error(three_and_a_half), EPSILON);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_par_reload_mesh)
{
  const std::string field_name = "test_field";
  // 1D mesh divided into 10 segments
  // mfem::Mesh mesh(10);
  // 2D mesh divided into triangles
  mfem::Mesh mesh(3, 3, 3, mfem::Element::TETRAHEDRON);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::ofstream fout("mesh_par_" + std::to_string(rank));
  parmesh.ParPrint(fout);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::ParFiniteElementSpace parfes(&parmesh, &fec);

  // The mesh must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(COLL_NAME, &parmesh, owns_mesh);

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

  std::ofstream fout_rel("mesh_par_reload_" + std::to_string(rank));
  dynamic_cast<mfem::ParMesh*>(sdc_reader.GetMesh())->ParPrint(fout_rel);

  int n_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
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
  }

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
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
