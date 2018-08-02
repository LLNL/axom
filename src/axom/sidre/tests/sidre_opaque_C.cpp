/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "gtest/gtest.h"

#include "axom/sidre/interface/sidre.h"

//------------------------------------------------------------------------------
// Some simple types and functions used in tests
//------------------------------------------------------------------------------

enum Centering { _Zone_,
                 _Node_,
                 _UnknownCentering_};

enum DType { _Double_,
             _Int_,
             _UnknownType_};

typedef struct
{
  int ilo;
  int ihi;
} AA_extent;

AA_extent* AA_extent_new(int lo, int hi)
{
  AA_extent* self = (AA_extent*) malloc( sizeof(AA_extent) );
  self->ilo = lo;
  self->ihi = hi;
  return self;
}

void AA_extent_delete(AA_extent* self)
{
  free(self);
}

int AA_get_num_pts(AA_extent* self, Centering cent)
{
  int retval = 0;

  switch ( cent )
  {
  case _Zone_:
    retval = (self->ihi - self->ilo + 1);
    break;
  case _Node_:
    retval = (self->ihi - self->ilo + 2);
    break;
  default:
    retval = -1;       // I know magic numbers are bad. So sue me.
  }

  return retval;
}

typedef struct
{
  Centering cent;
  DType type;
  int depth;
} AA_meshvar;


AA_meshvar* AA_meshvar_new(Centering cent, DType type, int depth)
{
  AA_meshvar* self = (AA_meshvar*) malloc( sizeof(AA_meshvar) );
  self->cent = cent;
  self->type = type;
  self->depth = depth;
  return self;
}

void AA_meshvar_delete(AA_meshvar* self)
{
  free(self);
}

int AA_get_num_vals(AA_meshvar* self, AA_extent* ext)
{
  return AA_get_num_pts(ext, self->cent) * self->depth;
}

//------------------------------------------------------------------------------
//
// Simple test that adds an opaque data object, retrieves it and checks if
// the retrieved object is in the expected state.
//
TEST(C_sidre_opaque,basic_inout)
{
  const int ihi_val = 9;

  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);

  SIDRE_group* problem_gp = SIDRE_group_create_group(root, "problem");

  AA_extent* ext = AA_extent_new(0, ihi_val);

  SIDRE_view* ext_view = SIDRE_group_create_view_external(problem_gp,
                                                          "ext", ext);

  bool test_external = SIDRE_view_is_external(ext_view);
  EXPECT_EQ(test_external, true);

  bool test_applied = SIDRE_view_is_applied(ext_view);
  EXPECT_EQ(test_applied, false);

  bool test_opaque = SIDRE_view_is_opaque(ext_view);
  EXPECT_EQ(test_opaque, true);

  AA_extent* test_extent = (AA_extent*) SIDRE_view_get_void_ptr(ext_view);
  int test_ihi = test_extent->ihi;

  EXPECT_EQ(test_ihi, ihi_val);

#if 1
  // Similar test with different view methods

  AA_extent* ext2 = AA_extent_new(0, 2 * ihi_val);

  SIDRE_view* ext2_view = SIDRE_group_create_view_empty(problem_gp,
                                                        "ext2");
  SIDRE_view_set_external_data_ptr_only(ext2_view, ext2);

  bool test_opaque2 = SIDRE_view_is_opaque(ext2_view);
  EXPECT_EQ(test_opaque2, true);

  AA_extent* test_extent2 =
    (AA_extent*) SIDRE_view_get_void_ptr(ext2_view);
  int test_ihi2 = test_extent2->ihi;

  EXPECT_EQ(test_ihi2, 2 * ihi_val);
#endif

  // clean up...
  AA_extent_delete(ext);
  SIDRE_datastore_delete(ds);
}

//------------------------------------------------------------------------------
//
// Test that adds "MeshVars" as opaque data objects, creates views for their
// data on each of two domains, allocates their data (based on centering,
// domain size, and depth), and then checks to if the allocated data
// lengths match the expected values.
//
TEST(C_sidre_opaque,meshvar)
{
  const int ilo_val[] = {0, 10};
  const int ihi_val[] = {9, 21};
  const std::string dom_name[] = { std::string("domain0"),
                                   std::string("domain1") };

  const int zone_var_depth = 1;
  const int node_var_depth = 2;

  SIDRE_datastore* ds = SIDRE_datastore_new();
  SIDRE_group* root = SIDRE_datastore_get_root(ds);

  SIDRE_group* problem_gp = SIDRE_group_create_group(root, "problem");

  // Add two different mesh vars to mesh var group
  SIDRE_group* meshvar_gp =
    SIDRE_group_create_group(problem_gp, "mesh_var");

  AA_meshvar* zone_mv = AA_meshvar_new(_Zone_, _Int_, zone_var_depth);
  SIDRE_view* zone_mv_view =
    SIDRE_group_create_view_external(meshvar_gp, "zone_mv", zone_mv);

  AA_meshvar* node_mv = AA_meshvar_new(_Node_, _Double_, node_var_depth);
  SIDRE_view* node_mv_view =
    SIDRE_group_create_view_external(meshvar_gp, "node_mv", node_mv);

  //
  // Create domain groups, add extents
  // Create data views for mesh var data on domains and allocate
  //
  for (int idom = 0 ; idom < 2 ; ++idom)
  {
    SIDRE_group* dom_gp =
      SIDRE_group_create_group(problem_gp, dom_name[idom].c_str());

    AA_extent* dom_ext = AA_extent_new(ilo_val[idom], ihi_val[idom]);
    SIDRE_group_create_view_external(dom_gp, "ext", dom_ext);

    AA_meshvar* zonemv =
      (AA_meshvar*) SIDRE_view_get_void_ptr(zone_mv_view);
    (void) SIDRE_group_create_view_and_allocate_nelems(dom_gp, "zone_data",
                                                       SIDRE_INT_ID,
                                                       AA_get_num_vals(
                                                         zonemv,
                                                         dom_ext));

    AA_meshvar* nodemv =
      (AA_meshvar*)  SIDRE_view_get_void_ptr(node_mv_view);
    (void) SIDRE_group_create_view_and_allocate_nelems(dom_gp, "node_data",
                                                       SIDRE_DOUBLE_ID, AA_get_num_vals(
                                                         nodemv,
                                                         dom_ext));

  }

//
//  Print datastore contents to see what's going on.
//
//  SIDRE_datastore_print(ds);


  //
  // Check data lengths
  //
  for (int idom = 0 ; idom < 2 ; ++idom)
  {
    SIDRE_group* dom_gp =
      SIDRE_group_get_group_from_name(problem_gp, dom_name[idom].c_str());
    SIDRE_view* ext_view =
      SIDRE_group_get_view_from_name(dom_gp, "ext");
    AA_extent* dom_ext = (AA_extent*) SIDRE_view_get_void_ptr(ext_view);

    AA_meshvar* zonemv = (AA_meshvar*) SIDRE_view_get_void_ptr(
      zone_mv_view);
    AA_meshvar* nodemv = (AA_meshvar*) SIDRE_view_get_void_ptr(
      node_mv_view);

    int num_zone_vals = AA_get_num_vals(zonemv, dom_ext);
    SIDRE_view* dom_zone_data_view = SIDRE_group_get_view_from_name(
      dom_gp, "zone_data");
    int test_num_zone_vals =
      SIDRE_view_get_num_elements(dom_zone_data_view);
    EXPECT_EQ(num_zone_vals, test_num_zone_vals);

    int num_node_vals = AA_get_num_vals(nodemv, dom_ext);
    SIDRE_view* dom_node_data_view = SIDRE_group_get_view_from_name(
      dom_gp, "node_data");
    int test_num_node_vals =
      SIDRE_view_get_num_elements(dom_node_data_view);
    EXPECT_EQ(num_node_vals, test_num_node_vals);

  }

  // clean up...
  AA_meshvar_delete(zone_mv);
  AA_meshvar_delete(node_mv);
  for (int idom = 0 ; idom < 2 ; ++idom)
  {
    SIDRE_group* dom_gp =
      SIDRE_group_get_group_from_name(problem_gp, dom_name[idom].c_str());
    SIDRE_view* ext_view =
      SIDRE_group_get_view_from_name(dom_gp, "ext");
    AA_extent* dom_ext = (AA_extent*) SIDRE_view_get_void_ptr(ext_view);
    AA_extent_delete(dom_ext);
  }
  SIDRE_datastore_delete(ds);
}
