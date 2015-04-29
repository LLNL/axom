#include "gtest/gtest.h"

#include "datastore/datastore.h"

using namespace sidre;
using namespace conduit;

//------------------------------------------------------------------------------
TEST(datastore_group,create_group)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataGroup *flds = root->createGroup("fields");
    
    DataBuffer *db0 = ds->createBuffer();
    DataBuffer *db1 = ds->createBuffer();
    
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
    
    flds->createView("u0",db0)->apply(DataType::uint64(10));
    flds->createView("f1",db1)->apply(DataType::float64(10));
    
    flds->getView("u0")->getNode().print_detailed();
    flds->getView("f1")->getNode().print_detailed();


    delete ds;
}

//------------------------------------------------------------------------------
TEST(datastore_group,group_name_collisions)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->getRoot()->createGroup("fields");
    flds->createViewAndBuffer("a");

    EXPECT_TRUE(flds->hasView("a"));

    delete ds;
}
//------------------------------------------------------------------------------
TEST(datastore_group,view_copy_move)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->getRoot()->createGroup("fields");

    flds->createViewAndBuffer("i0")->allocate(DataType::int32());
    flds->createViewAndBuffer("f0")->allocate(DataType::float32());
    flds->createViewAndBuffer("d0")->allocate(DataType::float64());

    (*flds->getView("i0")->getNode().as_int32_ptr())   = 1;
    (*flds->getView("f0")->getNode().as_float32_ptr()) = 100.0;
    (*flds->getView("d0")->getNode().as_float64_ptr()) = 3000.0;

    EXPECT_TRUE(flds->hasView("i0"));
    EXPECT_TRUE(flds->hasView("f0"));
    EXPECT_TRUE(flds->hasView("d0"));

    // test moving a view form feds7 to sub
    flds->createGroup("sub")->moveView(flds->getView("d0"));
    flds->print();
    EXPECT_FALSE(flds->hasView("d0"));
    EXPECT_TRUE(flds->hasGroup("sub"));
    EXPECT_TRUE(flds->getGroup("sub")->hasView("d0"));

    // check the data value
    float64 *d0_data =  flds->getGroup("sub")
                            ->getView("d0")
                            ->getNode().as_float64_ptr();
    EXPECT_NEAR(d0_data[0],3000.0,1e-12);
    
    // test copying a view from flds top sub
    flds->getGroup("sub")->copyView(flds->getView("i0"));

    flds->print();
    
    EXPECT_TRUE(flds->hasView("i0"));    
    EXPECT_TRUE(flds->getGroup("sub")->hasView("i0"));

    // we expect the actual data  pointers to be the same
    EXPECT_EQ(flds->getView("i0")->getNode().data_pointer(),
              flds->getGroup("sub")->getView("i0")->getNode().data_pointer());

    delete ds;
}

//------------------------------------------------------------------------------
TEST(datastore_group,groups_move_copy)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->getRoot()->createGroup("fields");

    DataGroup *ga = flds->createGroup("a");
    DataGroup *gb = flds->createGroup("b");
    DataGroup *gc = flds->createGroup("c");

    ga->createViewAndBuffer("i0")->allocate(DataType::int32());
    gb->createViewAndBuffer("f0")->allocate(DataType::float32());
    gc->createViewAndBuffer("d0")->allocate(DataType::float64());

    (*ga->getView("i0")->getNode().as_int32_ptr())   = 1;
    (*gb->getView("f0")->getNode().as_float32_ptr()) = 100.0;
    (*gc->getView("d0")->getNode().as_float64_ptr()) = 3000.0;

    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("b"));
    EXPECT_TRUE(flds->hasGroup("c"));

    //move "b" to a child of "sub"
    flds->createGroup("sub")->moveGroup(gb);

    flds->print();
    
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("sub"));
    EXPECT_TRUE(flds->hasGroup("c"));

    EXPECT_EQ(flds->getGroup("sub")->getGroup("b"),gb);

    delete ds;
}


TEST(datastore_group,create_destroy_view_and_buffer)
{
  DataStore * const ds = new DataStore();
  DataGroup * const grp = ds->getRoot()->createGroup("grp");

  std::string const viewName1 = "viewBuffer1";
  std::string const viewName2 = "viewBuffer2";

  DataView const * const view1 = grp->createViewAndBuffer(viewName1);
  DataView const * const view2 = grp->createViewAndBuffer(viewName2);

  EXPECT_TRUE(grp->hasView(viewName1));
  EXPECT_EQ( grp->getView(viewName1), view1 );

  EXPECT_TRUE(grp->hasView(viewName2));
  EXPECT_EQ( grp->getView(viewName2), view2 );

  IDType const bufferId1 = view1->getBuffer()->GetUID();

  grp->destroyViewAndBuffer(viewName1);


  EXPECT_FALSE(grp->hasView(viewName1));
  EXPECT_EQ(ds->getNumberOfBuffers(),1);

  DataBuffer const * const buffer1 = ds->getBuffer(bufferId1);
  bool buffValid = true;
  if( buffer1 == nullptr )
  {
    buffValid = false;
  }

  EXPECT_FALSE(buffValid);

  delete ds;
}



//------------------------------------------------------------------------------
TEST(datastore_group,save_restore_simple)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->getRoot()->createGroup("fields");

    DataGroup *ga = flds->createGroup("a");

    ga->createViewAndBuffer("i0")->allocate(DataType::int32());

    (*ga->getView("i0")->getNode().as_int32_ptr())   = 1;

    EXPECT_TRUE(ds->getRoot()->hasGroup("fields"));
    EXPECT_TRUE(ds->getRoot()->getGroup("fields")->hasGroup("a"));
    EXPECT_TRUE(ds->getRoot()->getGroup("fields")->getGroup("a")->hasView("i0"));
        
        
    ds->getRoot()->save("out_ds_group_save_restore_simple","conduit");
    
    ds->print();
    
    DataStore *ds2 = new DataStore();

    ds2->getRoot()->load("out_ds_group_save_restore_simple","conduit");
    
    flds = ds2->getRoot()->getGroup("fields");
    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_EQ(flds->getGroup("a")->getView("i0")->getNode().as_int32(),1);
    
    ds2->print();
    
    delete ds;    
    delete ds2;
    
}



//------------------------------------------------------------------------------
TEST(datastore_group,save_restore_complex)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->getRoot()->createGroup("fields");

    DataGroup *ga = flds->createGroup("a");
    DataGroup *gb = flds->createGroup("b");
    DataGroup *gc = flds->createGroup("c");

    ga->createViewAndBuffer("i0")->allocate(DataType::int32());
    gb->createViewAndBuffer("f0")->allocate(DataType::float32());
    gc->createViewAndBuffer("d0")->allocate(DataType::float64());

    (*ga->getView("i0")->getNode().as_int32_ptr())   = 1;
    (*gb->getView("f0")->getNode().as_float32_ptr()) = 100.0;
    (*gc->getView("d0")->getNode().as_float64_ptr()) = 3000.0;

    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("b"));
    EXPECT_TRUE(flds->hasGroup("c"));

    ds->getRoot()->save("out_ds_group_save_restore_complex","conduit");

    ds->print();

    DataStore *ds2 = new DataStore();


    ds2->getRoot()->load("out_ds_group_save_restore_complex","conduit");

    flds = ds2->getRoot()->getGroup("fields");
    // check that all sub groups exist
    EXPECT_TRUE(flds->hasGroup("a"));
    EXPECT_TRUE(flds->hasGroup("b"));
    EXPECT_TRUE(flds->hasGroup("c"));

    EXPECT_EQ(flds->getGroup("a")->getView("i0")->getNode().as_int32(),1);
    EXPECT_NEAR(flds->getGroup("b")->getView("f0")->getNode().as_float32(),100.0,  1e-12);
    EXPECT_NEAR(flds->getGroup("c")->getView("d0")->getNode().as_float64(),3000.0, 1e-12);

    ds2->print();

    delete ds;
    delete ds2;

}


