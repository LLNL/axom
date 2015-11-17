/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include "sidre/sidre.h"

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

AA_extent * AA_extent_new(int lo, int hi)
{
  AA_extent * self = (AA_extent *) malloc( sizeof(AA_extent) );
  self->ilo = lo;
  self->ihi = hi;
  return self;
}

void AA_extent_delete(AA_extent * self)
{
  free(self);
}

int AA_get_num_pts(AA_extent * self, Centering cent)
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


AA_meshvar * AA_meshvar_new(Centering cent, DType type, int depth)
{
  AA_meshvar * self = (AA_meshvar *) malloc( sizeof(AA_meshvar) );
  self->cent = cent;
  self->type = type;
  self->depth = depth;
  return self;
}

void AA_meshvar_delete(AA_meshvar * self)
{
  free(self);
}

int AA_get_num_vals(AA_meshvar * self, AA_extent * ext)
{
  return AA_get_num_pts(ext, self->cent) * self->depth;
}

//------------------------------------------------------------------------------
//
// Simple test that adds an opaque data object, retrieves it and checks if
// the retrieved object is in the expected state.
//
TEST(C_sidre_opaque,inout)
{
  const int ihi_val = 9;

  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);

  ATK_datagroup * problem_gp = ATK_datagroup_create_group(root, "problem");

  AA_extent * ext = AA_extent_new(0, ihi_val);

  ATK_dataview * ext_view = ATK_datagroup_create_opaque_view(problem_gp, "ext", ext);
#if 1
//  ATK_datagroup_CreateViewAndBuffer("ext");
//  ATK_datagroup_CreateOpaqueView("ext", ext);
//  ATK_datagroup_CreateView("ext", 0);
//  ATK_datagroup_MoveView(0);
//  ATK_datagroup_MoveView(ATK_datagroup_GetView(problem_gp, "ext"));
//  ATK_datagroup_CopyView(0);
//  ATK_datagroup_CopyView(ATK_datagroup_GetView(problem_gp, "ext"));
//  ATK_datagroup_AttachView(0);
//  ATK_datagroup_CopyView(ATK_datagroup_GetView(problem_gp, "ext"));
//  Can't do following: method is private...
//  DataView* v = ATK_datagroup_DetachView("ext");
//  std::cout << "view name = " << v->GetName() << std::endl;
//  ATK_datagroup_DestroyView("foo");
//  root->MoveGroup(problem_gp);
//  root->CopyGroup(problem_gp);
//  Can't do following: method is private...
//  root->DetachGroup("bar");
//  root->DestroyGroup("bar");
//  ATK_datagroup_get_view(2);
#endif

  bool test_opaque = ATK_dataview_is_opaque(ext_view);
  EXPECT_EQ(test_opaque, true);

  AA_extent * test_extent = (AA_extent *) ATK_dataview_get_opaque(ext_view);
  int test_ihi = test_extent->ihi;

  EXPECT_EQ(test_ihi, ihi_val);

  // clean up...
  AA_extent_delete(ext);
  ATK_datastore_delete(ds);
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

  ATK_datastore * ds = ATK_datastore_new();
  ATK_datagroup * root = ATK_datastore_get_root(ds);

  ATK_datagroup * problem_gp = ATK_datagroup_create_group(root, "problem");

  // Add two different mesh vars to mesh var group
  ATK_datagroup * meshvar_gp = ATK_datagroup_create_group(problem_gp, "mesh_var");
  AA_meshvar * zone_mv = AA_meshvar_new(_Zone_, _Int_, zone_var_depth);
  ATK_dataview * zone_mv_view = ATK_datagroup_create_opaque_view(meshvar_gp, "zone_mv", zone_mv);
  AA_meshvar * node_mv = AA_meshvar_new(_Node_, _Double_, node_var_depth);
  ATK_dataview * node_mv_view = ATK_datagroup_create_opaque_view(meshvar_gp, "node_mv", node_mv);

  //
  // Create domain groups, add extents
  // Create data views for mesh var data on domains and allocate
  //
  for (int idom = 0 ; idom < 2 ; ++idom)
  {
    ATK_datagroup * dom_gp = ATK_datagroup_create_group(problem_gp, dom_name[idom].c_str());
    AA_extent * dom_ext = AA_extent_new(ilo_val[idom], ihi_val[idom]);
    ATK_datagroup_create_opaque_view(dom_gp, "ext", dom_ext);

    AA_meshvar * zonemv = (AA_meshvar *) ATK_dataview_get_opaque(zone_mv_view);
    ATK_dataview * dom_zone_view = ATK_datagroup_create_view_and_buffer_simple(dom_gp, "zone_data");
    ATK_dataview_allocate_from_type(dom_zone_view, SIDRE_INT_ID, AA_get_num_vals(zonemv, dom_ext));

    AA_meshvar * nodemv = (AA_meshvar *)  ATK_dataview_get_opaque(node_mv_view);
    ATK_dataview * dom_node_view = ATK_datagroup_create_view_and_buffer_simple(dom_gp, "node_data");
    ATK_dataview_allocate_from_type(dom_node_view, SIDRE_DOUBLE_ID, AA_get_num_vals(nodemv, dom_ext));

  }

//
//  Print datastore contents to see what's going on.
//
//  ATK_datastore_print(ds);


  //
  // Check data lengths
  //
  for (int idom = 0 ; idom < 2 ; ++idom)
  {
    ATK_datagroup * dom_gp = ATK_datagroup_get_group(problem_gp, dom_name[idom].c_str());
    ATK_dataview * ext_view = ATK_datagroup_get_view_from_name(dom_gp, "ext");
    AA_extent * dom_ext = (AA_extent *) ATK_dataview_get_opaque(ext_view);

    AA_meshvar * zonemv = (AA_meshvar *) ATK_dataview_get_opaque(zone_mv_view);
    AA_meshvar * nodemv = (AA_meshvar *) ATK_dataview_get_opaque(node_mv_view);

    int num_zone_vals = AA_get_num_vals(zonemv, dom_ext);
    ATK_dataview * dom_zone_data_view = ATK_datagroup_get_view_from_name(dom_gp, "zone_data");
    int test_num_zone_vals = ATK_dataview_get_number_of_elements(dom_zone_data_view);
    EXPECT_EQ(num_zone_vals, test_num_zone_vals);

    int num_node_vals = AA_get_num_vals(nodemv, dom_ext);
    ATK_dataview * dom_node_data_view = ATK_datagroup_get_view_from_name(dom_gp, "node_data");
    int test_num_node_vals = ATK_dataview_get_number_of_elements(dom_node_data_view);
    EXPECT_EQ(num_node_vals, test_num_node_vals);

  }

  // clean up...
  AA_meshvar_delete(zone_mv);
  AA_meshvar_delete(node_mv);
  for (int idom = 0 ; idom < 2 ; ++idom)
  {
    ATK_datagroup * dom_gp = ATK_datagroup_get_group(problem_gp, dom_name[idom].c_str());
    ATK_dataview * ext_view = ATK_datagroup_get_view_from_name(dom_gp, "ext");
    AA_extent * dom_ext = (AA_extent *) ATK_dataview_get_opaque(ext_view);
    AA_extent_delete(dom_ext);
  }
  ATK_datastore_delete(ds);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;   // create & initialize test logger,
  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
