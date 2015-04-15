#include "gtest/gtest.h"

#include "datastore/datastore.h"

using namespace DataStoreNS;
using namespace conduit;

//------------------------------------------------------------------------------
TEST(datastore_group,create_group)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->GetRoot();
    DataGroup *flds = root->CreateGroup("fields");
    
    DataBuffer *db0 = ds->CreateBuffer();
    DataBuffer *db1 = ds->CreateBuffer();
    
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


    delete ds;
}

//------------------------------------------------------------------------------
TEST(datastore_group,group_name_collisions)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->GetRoot()->CreateGroup("fields");
    flds->CreateViewAndBuffer("a");

    EXPECT_TRUE(flds->HasChild("a"));

    delete ds;
}
//------------------------------------------------------------------------------
TEST(datastore_group,view_copy_move)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->GetRoot()->CreateGroup("fields");

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

    delete ds;
}

//------------------------------------------------------------------------------
TEST(datastore_group,groups_move_copy)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->GetRoot()->CreateGroup("fields");

    DataGroup *ga = flds->CreateGroup("a");
    DataGroup *gb = flds->CreateGroup("b");
    DataGroup *gc = flds->CreateGroup("c");

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

    delete ds;
}


TEST(datastore_group,create_destroy_view_and_buffer)
{
  DataStore * const ds = new DataStore();
  DataGroup * const grp = ds->GetRoot()->CreateGroup("grp");

  std::string const viewName1 = "viewBuffer1";
  std::string const viewName2 = "viewBuffer2";

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

  DataBuffer const * const buffer1 = ds->GetBuffer(bufferId1);
  bool buffValid = true;
  if( buffer1 == nullptr )
  {
    buffValid = false;
  }

  EXPECT_FALSE(buffValid);

  delete ds;
}

