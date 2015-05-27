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
TEST(C_sidre_group,create_group)
{
    DS_datastore *ds = DS_datastore_new();
    DS_datagroup *root = DS_datastore_get_root(ds);
    DS_datagroup *flds = DS_datagroup_create_group(root, "fields");
    EXPECT_TRUE(flds != NULL);

#if 0
    DS_databuffer *db0 = DS_datastore_create_buffer(ds);
    DS_databuffer *db1 = DS_datastore_create_buffer(ds);
    
    db0->Declare(DataType::uint64(10));
    db0->Allocate();
    uint64 *db0_ptr = db0->GetNode().as_uint64_ptr();
    
    db1->Declare(DataType::float64(10));
    db1->Allocate();
    float64 *db1_ptr = db1->GetNode().as_float64_ptr();
    
    for(int i=0;i<10;i++)
    {
        db0_ptr[i] = i;
        db1_ptr[i] = i*i;
    }
    
    flds->CreateView("u0",db0)->Apply(DataType::uint64(10));
    flds->CreateView("f1",db1)->Apply(DataType::float64(10));
    
    flds->GetView("u0")->GetNode().print_detailed();
    flds->GetView("f1")->GetNode().print_detailed();

#endif
    DS_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,group_name_collisions)
{
    DS_datastore *ds = DS_datastore_new();
    DS_datagroup *root = DS_datastore_get_root(ds);
    DS_datagroup *flds = DS_datagroup_create_group(root, "fields");

    DS_datagroup_create_view_and_buffer(flds, "a");

    EXPECT_TRUE(flds != NULL);
    //    EXPECT_TRUE(flds->HasView("a"));

    DS_datastore_delete(ds);
}
//------------------------------------------------------------------------------
TEST(C_sidre_group,view_copy_move)
{
    DS_datastore *ds = DS_datastore_new();
    DS_datagroup *root = DS_datastore_get_root(ds);
    DS_datagroup *flds = DS_datagroup_create_group(root, "fields");

    EXPECT_TRUE(flds != NULL);
#if 0
    flds->CreateViewAndBuffer("i0")->Allocate(DataType::int32());
    flds->CreateViewAndBuffer("f0")->Allocate(DataType::float32());
    flds->CreateViewAndBuffer("d0")->Allocate(DataType::float64());

    (*flds->GetView("i0")->GetNode().as_int32_ptr())   = 1;
    (*flds->GetView("f0")->GetNode().as_float32_ptr()) = 100.0;
    (*flds->GetView("d0")->GetNode().as_float64_ptr()) = 3000.0;

    EXPECT_TRUE(flds->HasView("i0"));
    EXPECT_TRUE(flds->HasView("f0"));
    EXPECT_TRUE(flds->HasView("d0"));

    // test moving a view form feds7 to sub
    flds->CreateGroup("sub")->MoveView(flds->GetView("d0"));
    flds->Print();
    EXPECT_FALSE(flds->HasView("d0"));
    EXPECT_TRUE(flds->HasGroup("sub"));
    EXPECT_TRUE(flds->GetGroup("sub")->HasView("d0"));

    // check the data value
    float64 *d0_data =  flds->GetGroup("sub")
                            ->GetView("d0")
                            ->GetNode().as_float64_ptr();
    EXPECT_NEAR(d0_data[0],3000.0,1e-12);
    
    // test copying a view from flds top sub
    flds->GetGroup("sub")->CopyView(flds->GetView("i0"));

    flds->Print();
    
    EXPECT_TRUE(flds->HasView("i0"));    
    EXPECT_TRUE(flds->GetGroup("sub")->HasView("i0"));

    // we expect the actual data  pointers to be the same
    EXPECT_EQ(flds->GetView("i0")->GetNode().data_pointer(),
              flds->GetGroup("sub")->GetView("i0")->GetNode().data_pointer());
#endif

    DS_datastore_delete(ds);
}

//------------------------------------------------------------------------------
TEST(C_sidre_group,groups_move_copy)
{
    DS_datastore *ds = DS_datastore_new();
    DS_datagroup *root = DS_datastore_get_root(ds);
    DS_datagroup *flds = DS_datagroup_create_group(root, "fields");
    EXPECT_TRUE(flds != NULL);

#if 0
    DS_datagroup *ga = DS_datagroup_create_group(flds, "a");
    DS_datagroup *gb = DS_datagroup_create_group(flds, "b");
    DS_datagroup *gc = DS_datagroup_create_group(flds, "c");

    ga->CreateViewAndBuffer("i0")->Allocate(DataType::int32());
    gb->CreateViewAndBuffer("f0")->Allocate(DataType::float32());
    gc->CreateViewAndBuffer("d0")->Allocate(DataType::float64());

    (*ga->GetView("i0")->GetNode().as_int32_ptr())   = 1;
    (*gb->GetView("f0")->GetNode().as_float32_ptr()) = 100.0;
    (*gc->GetView("d0")->GetNode().as_float64_ptr()) = 3000.0;

    // check that all sub groups exist
    EXPECT_TRUE(flds->HasGroup("a"));
    EXPECT_TRUE(flds->HasGroup("b"));
    EXPECT_TRUE(flds->HasGroup("c"));

    //move "b" to a child of "sub"
    flds->CreateGroup("sub")->MoveGroup(gb);

    flds->Print();
    
    EXPECT_TRUE(flds->HasGroup("a"));
    EXPECT_TRUE(flds->HasGroup("sub"));
    EXPECT_TRUE(flds->HasGroup("c"));

    EXPECT_EQ(flds->GetGroup("sub")->GetGroup("b"),gb);
#endif

    DS_datastore_delete(ds);
}


TEST(C_sidre_group,create_destroy_view_and_buffer)
{
  DS_datastore * const ds = DS_datastore_new();
  DS_datagroup *root = DS_datastore_get_root(ds);
  DS_datagroup *grp = DS_datagroup_create_group(root, "grp");
  EXPECT_TRUE(grp != NULL);

  //char *viewName1 = "viewBuffer1";
  //char *viewName2 = "viewBuffer2";

#if 0
  DataView const * const view1 = grp->CreateViewAndBuffer(viewName1);
  DataView const * const view2 = grp->CreateViewAndBuffer(viewName2);

  EXPECT_TRUE(grp->HasView(viewName1));
  EXPECT_EQ( grp->GetView(viewName1), view1 );

  EXPECT_TRUE(grp->HasView(viewName2));
  EXPECT_EQ( grp->GetView(viewName2), view2 );

  IDType const bufferId1 = view1->GetBuffer()->GetUID();

  grp->DestroyViewAndBuffer(viewName1);


  EXPECT_FALSE(grp->HasView(viewName1));
  EXPECT_EQ(ds->GetNumberOfBuffers(),1);

  DS_databuffer const * const buffer1 = ds->GetBuffer(bufferId1);
  bool buffValid = true;
  if( buffer1 == nullptr )
  {
    buffValid = false;
  }

  EXPECT_FALSE(buffValid);
#endif

  DS_datastore_delete(ds);
}



//------------------------------------------------------------------------------
TEST(C_sidre_group,save_restore_simple)
{
    DS_datastore *ds = DS_datastore_new();
    DS_datagroup *root = DS_datastore_get_root(ds);
    DS_datagroup *flds = DS_datagroup_create_group(root, "fields");
    EXPECT_TRUE(flds != NULL);

#if 0
    DS_datagroup *ga = DS_datagroup_create_group(flds, "a");

    ga->CreateViewAndBuffer("i0")->Allocate(DataType::int32());

    (*ga->GetView("i0")->GetNode().as_int32_ptr())   = 1;

    EXPECT_TRUE(ds->GetRoot()->HasGroup("fields"));
    EXPECT_TRUE(ds->GetRoot()->GetGroup("fields")->HasGroup("a"));
    EXPECT_TRUE(ds->GetRoot()->GetGroup("fields")->GetGroup("a")->HasView("i0"));
        
        
    ds->GetRoot()->save("out_sidre_group_save_restore_simple","conduit");
    
    ds->Print();
    
    DS_datastore *ds2 = DS_datastore_new();

    ds2->GetRoot()->load("out_sidre_group_save_restore_simple","conduit");
    
    flds = ds2->GetRoot()->GetGroup("fields");
    // check that all sub groups exist
    EXPECT_TRUE(flds->HasGroup("a"));
    EXPECT_EQ(flds->GetGroup("a")->GetView("i0")->GetNode().as_int32(),1);
    
    ds2->Print();

    DS_datastore_delete(ds2);
#endif    
    DS_datastore_delete(ds);
    
}



//------------------------------------------------------------------------------
TEST(C_sidre_group,save_restore_complex)
{
    DS_datastore *ds = DS_datastore_new();
    DS_datagroup *root = DS_datastore_get_root(ds);
    DS_datagroup *flds = DS_datagroup_create_group(root, "fields");
    EXPECT_TRUE(flds != NULL);

#if 0
    DS_datagroup *ga = DS_datagroup_create_group(flds, "a");
    DS_datagroup *gb = DS_datagroup_create_group(flds, "b");
    DS_datagroup *gc = DS_datagroup_create_group(flds, "c");


    ga->CreateViewAndBuffer("i0")->Allocate(DataType::int32());
    gb->CreateViewAndBuffer("f0")->Allocate(DataType::float32());
    gc->CreateViewAndBuffer("d0")->Allocate(DataType::float64());

    (*ga->GetView("i0")->GetNode().as_int32_ptr())   = 1;
    (*gb->GetView("f0")->GetNode().as_float32_ptr()) = 100.0;
    (*gc->GetView("d0")->GetNode().as_float64_ptr()) = 3000.0;

    // check that all sub groups exist
    EXPECT_TRUE(flds->HasGroup("a"));
    EXPECT_TRUE(flds->HasGroup("b"));
    EXPECT_TRUE(flds->HasGroup("c"));

    ds->GetRoot()->save("out_sidre_group_save_restore_complex","conduit");

    ds->Print();

    DS_datastore *ds2 = DS_datastore_new();


    ds2->GetRoot()->load("out_sidre_group_save_restore_complex","conduit");

    flds = ds2->GetRoot()->GetGroup("fields");
    // check that all sub groups exist
    EXPECT_TRUE(flds->HasGroup("a"));
    EXPECT_TRUE(flds->HasGroup("b"));
    EXPECT_TRUE(flds->HasGroup("c"));

    EXPECT_EQ(flds->GetGroup("a")->GetView("i0")->GetNode().as_int32(),1);
    EXPECT_NEAR(flds->GetGroup("b")->GetView("f0")->GetNode().as_float32(),100.0,  1e-12);
    EXPECT_NEAR(flds->GetGroup("c")->GetView("d0")->GetNode().as_float64(),3000.0, 1e-12);

    ds2->Print();
    DS_datastore_delete(ds2);
#endif

    DS_datastore_delete(ds);

}

