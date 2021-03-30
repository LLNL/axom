// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#ifndef AXOM_USE_MFEM
  #error This file requires MFEM
#endif

#include "mfem.hpp"

#include "gtest/gtest.h"

#include "axom/sidre/core/sidre.hpp"
#include "axom/sidre/core/MFEMSidreDataCollection.hpp"

using axom::sidre::Group;
using axom::sidre::MFEMSidreDataCollection;

const double EPSILON = 1.0e-6;

std::string testName()
{
  return ::testing::UnitTest::GetInstance()->current_test_info()->name();
}

TEST(sidre_datacollection, dc_alloc_no_mesh)
{
  MFEMSidreDataCollection sdc(testName());
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_owning_mesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc(testName(), &mesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_nonowning_mesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  bool owns_mesh = false;
  MFEMSidreDataCollection sdc(testName(), &mesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_register_empty_field)
{
  MFEMSidreDataCollection sdc(testName());
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

  MFEMSidreDataCollection sdc(testName(), &mesh);
  sdc.RegisterField("test_field", &gf);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_update_state)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  MFEMSidreDataCollection sdc(testName(), &mesh);

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
  MFEMSidreDataCollection sdc(testName(), &mesh);

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
  MFEMSidreDataCollection sdc_writer(testName(), &mesh, owns_mesh);
  mfem::GridFunction gf_write(&fes, nullptr);

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, &gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write.ProjectCoefficient(three_and_a_half);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
#endif

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  // No mesh is used here
  MFEMSidreDataCollection sdc_reader(testName());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
#endif

  sdc_reader.Load();

  // No need to reregister, it already exists
  auto gf_read = sdc_reader.GetField(field_name);

  // Make sure the gridfunction was actually read in
  EXPECT_LT(gf_read->ComputeL2Error(three_and_a_half), EPSILON);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_reload_gf_vdim)
{
  const std::string field_name = "test_field";
  const int vdim = 2;
  // 2D mesh divided into triangles
  mfem::Mesh mesh(10, 10, mfem::Element::TRIANGLE);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::FiniteElementSpace fes(&mesh, &fec, vdim, mfem::Ordering::byVDIM);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(testName(), &mesh, owns_mesh);
  mfem::GridFunction gf_write(&fes, nullptr);

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, &gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write.ProjectCoefficient(three_and_a_half);

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
#endif

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  // No mesh is used here
  MFEMSidreDataCollection sdc_reader(testName());

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
#endif

  sdc_reader.Load();

  // No need to reregister, it already exists
  auto gf_read = sdc_reader.GetField(field_name);

  // Make sure the gridfunction was actually read in
  EXPECT_LT(gf_read->ComputeL2Error(three_and_a_half), EPSILON);
  EXPECT_EQ(gf_read->FESpace()->GetOrdering(), mfem::Ordering::byVDIM);
  EXPECT_EQ(gf_read->FESpace()->GetVDim(), vdim);

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
  MFEMSidreDataCollection sdc_writer(testName(), &mesh, owns_mesh);
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_writer.SetComm(MPI_COMM_WORLD);
#endif

  EXPECT_TRUE(sdc_writer.verifyMeshBlueprint());

  // Save some basic info about the mesh
  const int n_verts = sdc_writer.GetMesh()->GetNV();
  const int n_ele = sdc_writer.GetMesh()->GetNE();
  const int n_bdr_ele = sdc_writer.GetMesh()->GetNBE();

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  // No mesh is used here to construct as it will be read in
  MFEMSidreDataCollection sdc_reader(testName());
#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  sdc_reader.SetComm(MPI_COMM_WORLD);
#endif

  sdc_reader.Load();

  // Make sure the mesh was actually reconstructed
  EXPECT_EQ(sdc_reader.GetMesh()->GetNV(), n_verts);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNE(), n_ele);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNBE(), n_bdr_ele);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_qf_reload)
{
  //Set up a small mesh and a couple of grid function on that mesh
  mfem::Mesh* mesh =
    new mfem::Mesh(2, 3, mfem::Element::QUADRILATERAL, 0, 2.0, 3.0);
  mfem::FiniteElementCollection* fec = new mfem::LinearFECollection;
  mfem::FiniteElementSpace* fespace = new mfem::FiniteElementSpace(mesh, fec);

  const int intOrder = 3;
  const int qs_vdim = 1;
  const int qv_vdim = 2;

  mfem::QuadratureSpace* qspace = new mfem::QuadratureSpace(mesh, intOrder);
  // We want Sidre to allocate the data for us and to own the internal data.
  // If we don't do the below then the data collection doesn't save off the data
  // structure.
  mfem::QuadratureFunction* qs =
    new mfem::QuadratureFunction(qspace, nullptr, qs_vdim);
  mfem::QuadratureFunction* qv =
    new mfem::QuadratureFunction(qspace, nullptr, qv_vdim);

  int Nq = qs->Size();

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(COLL_NAME, mesh, owns_mesh);

  // sdc_writer owns the quadrature function fields and the underlying data
  sdc_writer.RegisterQField("qs", qs);
  sdc_writer.RegisterQField("qv", qv);

  // The data needs to be instantiated before we save it off
  for(int i = 0; i < Nq; ++i)
  {
    (*qs)(i) = double(i);
    (*qv)(2 * i + 0) = double(i);
    (*qv)(2 * i + 1) = double(Nq - i - 1);
  }

  sdc_writer.SetPrefixPath("/tmp/dc_qf_test");
  sdc_writer.SetCycle(5);
  sdc_writer.SetTime(8.0);
  sdc_writer.Save();

  MFEMSidreDataCollection sdc_reader(COLL_NAME);
  sdc_reader.SetPrefixPath("/tmp/dc_qf_test");
  sdc_reader.Load(sdc_writer.GetCycle());

  const int order_qs = sdc_reader.GetQFieldOrder("qs");
  const int vdim_qs = sdc_reader.GetQFieldVDim("qs");
  double* data_qs = sdc_reader.GetQFieldData("qs");

  const int order_qv = sdc_reader.GetQFieldOrder("qv");
  const int vdim_qv = sdc_reader.GetQFieldVDim("qv");

  // order_qs should also equal order_qv in this trivial case
  EXPECT_TRUE(order_qs == intOrder);
  EXPECT_TRUE(order_qv == intOrder);

  EXPECT_TRUE(vdim_qs == qs_vdim);
  EXPECT_TRUE(vdim_qv == qv_vdim);

  mfem::QuadratureSpace qspace_new(mesh, order_qs);

  // Just showing different ways in which you could restore a QuadratureFunction
  mfem::QuadratureFunction qs_new(&qspace_new, data_qs, vdim_qs);
  mfem::QuadratureFunction qv_new(&qspace_new,
                                  sdc_reader.GetQFieldData("qv"),
                                  sdc_reader.GetQFieldVDim("qv"));

  qs_new -= *qs;
  qv_new -= *qv;

  EXPECT_TRUE(qs_new.Norml2() < 1e-15);
  EXPECT_TRUE(qv_new.Norml2() < 1e-15);

  delete fespace;
  delete qs;
  delete qv;
  delete qspace;
  delete mesh;
}

#if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)

TEST(sidre_datacollection, dc_alloc_owning_parmesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  bool owns_mesh = true;
  MFEMSidreDataCollection sdc(testName(), &parmesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_alloc_nonowning_parmesh)
{
  // 1D mesh divided into 10 segments
  mfem::Mesh mesh(10);
  mfem::ParMesh parmesh(MPI_COMM_WORLD, mesh);
  bool owns_mesh = false;
  MFEMSidreDataCollection sdc(testName(), &parmesh, owns_mesh);
  EXPECT_TRUE(sdc.verifyMeshBlueprint());
}

struct ParMeshGroupData
{
  int n_verts = 0;
  int n_edges = 0;
  int n_triangles = 0;
  int n_quads = 0;
  bool operator==(const ParMeshGroupData& other) const
  {
    return (n_verts == other.n_verts) && (n_edges == other.n_edges) &&
      (n_triangles == other.n_triangles) && (n_quads == other.n_quads);
  }
};

static std::vector<ParMeshGroupData> getGroupData(const mfem::ParMesh& parmesh)
{
  const int dim = parmesh.Dimension();
  // Zeroth group doesn't matter
  std::vector<ParMeshGroupData> result(parmesh.GetNGroups() - 1);
  // remove when marked const in MFEM
  auto& non_const_parmesh = const_cast<mfem::ParMesh&>(parmesh);

  for(std::size_t i = 1; i <= result.size(); i++)
  {
    result[i - 1].n_verts = non_const_parmesh.GroupNVertices(i);
    if(dim >= 2)
    {
      result[i - 1].n_edges = non_const_parmesh.GroupNEdges(i);
      if(dim >= 3)
      {
        result[i - 1].n_triangles = non_const_parmesh.GroupNTriangles(i);
        result[i - 1].n_quads = non_const_parmesh.GroupNQuadrilaterals(i);
      }
    }
  }

  return result;
}

/**
 * @brief Helper method for testing that a parallel mesh is reconstructed correctly
 * @param [in] base_mesh The serial mesh object to distribute, save, and then reload
 * @param [in] part_method The partitioning method to use - in [0, 5]
 */
static void testParallelMeshReload(mfem::Mesh& base_mesh,
                                   const int part_method = 1)
{
  mfem::ParMesh parmesh(MPI_COMM_WORLD, base_mesh, nullptr, part_method);

  mfem::H1_FECollection fec(1, base_mesh.Dimension());
  mfem::ParFiniteElementSpace parfes(&parmesh, &fec);

  // The mesh must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  const bool owns_mesh = true;
  MFEMSidreDataCollection sdc_writer(testName(), &parmesh, owns_mesh);

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

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  MFEMSidreDataCollection sdc_reader(testName());

  // Needs to be set "manually" in order for everything to be loaded in properly
  sdc_reader.SetComm(MPI_COMM_WORLD);
  sdc_reader.Load();

  // Make sure the mesh was actually reconstructed
  EXPECT_EQ(sdc_reader.GetMesh()->GetNV(), n_verts);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNE(), n_ele);
  EXPECT_EQ(sdc_reader.GetMesh()->GetNBE(), n_bdr_ele);

  // Make sure the ParMesh was reconstructed - even on one rank, this
  // will still be a ParMesh
  auto reader_pmesh = dynamic_cast<mfem::ParMesh*>(sdc_reader.GetMesh());
  ASSERT_NE(reader_pmesh, nullptr);
  EXPECT_EQ(reader_pmesh->GetNGroups(), n_groups);
  EXPECT_EQ(reader_pmesh->GetNFaceNeighbors(), n_face_neighbors);
  EXPECT_EQ(reader_pmesh->GetNSharedFaces(), n_shared_faces);

  auto reader_group_data = getGroupData(*reader_pmesh);
  EXPECT_EQ(writer_group_data, reader_group_data);
  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

/**
 * @brief Wrapper for testing mesh reconstruction across all partitioning methods
 * @param [in] base_mesh The serial mesh object to distribute, save, and then reload
 */
static void testParallelMeshReloadAllPartitionings(mfem::Mesh& base_mesh)
{
  static constexpr int MAX_PART_METHOD =
    5;  // MFEM supports partition methods [0, 5]
  for(int part_method = 0; part_method <= MAX_PART_METHOD; part_method++)
  {
    testParallelMeshReload(base_mesh, part_method);
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
  MFEMSidreDataCollection sdc_writer(testName(), &parmesh, owns_mesh);

  // The mesh and field(s) must be owned by Sidre to properly manage data in case of
  // a simulated restart (save -> load)
  mfem::ParGridFunction gf_write(&parfes, static_cast<double*>(nullptr));

  // Register to allocate storage internally, then write to it
  sdc_writer.RegisterField(field_name, &gf_write);

  mfem::ConstantCoefficient three_and_a_half(3.5);
  gf_write.ProjectCoefficient(three_and_a_half);

  sdc_writer.SetCycle(0);
  sdc_writer.Save();

  MFEMSidreDataCollection sdc_reader(testName());

  // Needs to be set "manually" in order for everything to be loaded in properly
  sdc_reader.SetComm(MPI_COMM_WORLD);
  sdc_reader.Load();

  auto gf_read = sdc_reader.GetField(field_name);
  EXPECT_TRUE(dynamic_cast<mfem::ParGridFunction*>(gf_read));

  // Make sure the gridfunction was actually read in
  EXPECT_LT(gf_read->ComputeL2Error(three_and_a_half), EPSILON);

  EXPECT_TRUE(sdc_reader.verifyMeshBlueprint());
}

TEST(sidre_datacollection, dc_par_reload_mesh_1D_small)
{
  // 1D mesh divided into segments
  mfem::Mesh mesh(10);
  testParallelMeshReloadAllPartitionings(mesh);
}

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
  mfem::Mesh mesh(3, 3, 3, mfem::Element::HEXAHEDRON);
  testParallelMeshReloadAllPartitionings(mesh);
}

TEST(sidre_datacollection, dc_par_reload_mesh_3D_medium_hex)
{
  // 3D mesh divided into hexahedra
  mfem::Mesh mesh(10, 10, 10, mfem::Element::HEXAHEDRON);
  testParallelMeshReloadAllPartitionings(mesh);
}

  //----------------------------------------------------------------------
  #include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}

#endif  // defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
