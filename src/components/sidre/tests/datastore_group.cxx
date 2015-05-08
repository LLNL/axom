#include "gtest/gtest.h"

#include "sidre/sidre.hpp"

using asctoolkit::sidre::DataBuffer;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataView;
using asctoolkit::common::IDType;

using namespace conduit;

// API coverage tests
// Each test should be documented with the interface functions being tested

//------------------------------------------------------------------------------
// getName()
//------------------------------------------------------------------------------
TEST(datastore_group,get_name)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataGroup *group = root->createGroup("test");
 
    EXPECT_TRUE(group->getName() == std::string("test") );

    delete ds;
}

//------------------------------------------------------------------------------
// getParent()
//------------------------------------------------------------------------------
TEST(datastore_group,get_parent)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataGroup *parent = root->createGroup("parent");
    DataGroup *child = parent->createGroup("child");
 
    EXPECT_TRUE( child->getParent() == parent );

    delete ds;
}

//------------------------------------------------------------------------------
// Verify getDatastore()
//------------------------------------------------------------------------------
TEST(datastore_group,get_datastore)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataGroup *group = root->createGroup("parent");
 
    EXPECT_TRUE( group->getDataStore() == ds );

    DataStore const * const_ds = group->getDataStore();
    EXPECT_TRUE( const_ds == ds );

    delete ds;
}

//------------------------------------------------------------------------------
// hasGroup()
//------------------------------------------------------------------------------
TEST(datastore_group,has_child)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
 
    DataGroup *parent = root->createGroup("parent");
    DataGroup *child = parent->createGroup("child");
    EXPECT_TRUE( child->getParent() == parent );

    EXPECT_TRUE( parent->hasGroup("child") );

    delete ds;
}

//------------------------------------------------------------------------------
// createViewAndBuffer()
// destroyViewAndBuffer()
// hasView()
//------------------------------------------------------------------------------
TEST(datastore_group,create_destroy_has_viewbuffer)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataGroup *group = root->createGroup("parent");

    DataView *view = group->createViewAndBuffer("view");
    EXPECT_TRUE( group->getParent() == root );
    EXPECT_TRUE( view->hasBuffer() );
 
    EXPECT_TRUE( group->hasView("view") );

    group->destroyViewAndBuffer("view");

    EXPECT_FALSE( group->hasView("view") );

    delete ds;
}

//------------------------------------------------------------------------------
// createGroup()
// destroyGroup()
// hasGroup()
//------------------------------------------------------------------------------
TEST(datastore_group,create_destroy_has_group)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataGroup *group = root->createGroup("group");
    EXPECT_TRUE( group->getParent() == root );
    
    EXPECT_TRUE( root->hasGroup("group") );


    root->destroyGroup("group");
    EXPECT_FALSE( root->hasGroup("group") );

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

  IDType const bufferId1 = view1->getBuffer()->getUID();

  grp->destroyViewAndBuffer(viewName1);


  EXPECT_FALSE(grp->hasView(viewName1));
  EXPECT_EQ(ds->getNumberOfBuffers(),1u);

  DataBuffer const * const buffer1 = ds->getBuffer(bufferId1);
  bool buffValid = true;
  if( buffer1 == ATK_NULLPTR )
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


