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
TEST(datastore_group,group_name_collisons)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->GetRoot()->CreateGroup("fields");
    flds->CreateView("a");
    ASSERT_THROW(flds->CreateView("a"),std::exception);
    

    delete ds;
}

//------------------------------------------------------------------------------
TEST(datastore_group,attach_detach_views)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->GetRoot()->CreateGroup("fields");
    
    flds->CreateView("i0")->Allocate(DataType::int32());
    flds->CreateView("f0")->Allocate(DataType::float32());
    flds->CreateView("d0")->Allocate(DataType::float32());
    
    (*flds->GetView("i0")->GetNode().as_int32_ptr())   = 1;
    (*flds->GetView("f0")->GetNode().as_float32_ptr()) = 100.0;
    (*flds->GetView("d0")->GetNode().as_float64_ptr()) = 3000.0;
    
    EXPECT_TRUE(flds->HasView("i0"));
    EXPECT_TRUE(flds->HasView("f0"));
    EXPECT_TRUE(flds->HasView("d0"));
    
    DataView *f0 = flds->DetachView("f0");

    EXPECT_TRUE(flds->HasView("i0"));
    EXPECT_FALSE(flds->HasView("f0"));
    EXPECT_TRUE(flds->HasView("d0"));

    flds->AttachView(f0);
        
    EXPECT_TRUE(flds->HasView("i0"));
    EXPECT_TRUE(flds->HasView("f0"));
    EXPECT_TRUE(flds->HasView("d0"));
    
    flds->CreateGroup("sub")->AttachView(flds->DetachView("d0"));

    flds->Print();
    EXPECT_FALSE(flds->HasView("d0"));
    EXPECT_TRUE(flds->HasGroup("sub"));
    EXPECT_TRUE(flds->GetGroup("sub")->HasView("d0"));


    delete ds;
}


//------------------------------------------------------------------------------
TEST(datastore_group,attach_detach_groups)
{
    DataStore *ds = new DataStore();
    DataGroup *flds = ds->GetRoot()->CreateGroup("fields");
    
    DataGroup *ga = flds->CreateGroup("a");
    DataGroup *gb = flds->CreateGroup("b");
    DataGroup *gc = flds->CreateGroup("c");
    
    ga->CreateView("i0")->Allocate(DataType::int32());
    gb->CreateView("f0")->Allocate(DataType::float32());
    gc->CreateView("d0")->Allocate(DataType::float32());
    
    (*ga->GetView("i0")->GetNode().as_int32_ptr())   = 1;
    (*gb->GetView("f0")->GetNode().as_float32_ptr()) = 100.0;
    (*gc->GetView("d0")->GetNode().as_float64_ptr()) = 3000.0;

    // check that all sub groups exist
    EXPECT_TRUE(flds->HasGroup("a"));
    EXPECT_TRUE(flds->HasGroup("b"));
    EXPECT_TRUE(flds->HasGroup("c"));
    
    // detach b
    flds->DetachGroup("b");
    
    // make sure b isn't a sub of flds
    EXPECT_TRUE(flds->HasGroup("a"));
    EXPECT_FALSE(flds->HasGroup("b"));
    EXPECT_TRUE(flds->HasGroup("c"));
    
    // reattach b
    flds->AttachGroup(gb);
    
    // check that all sub groups exist
    EXPECT_TRUE(flds->HasGroup("a"));
    EXPECT_TRUE(flds->HasGroup("b"));
    EXPECT_TRUE(flds->HasGroup("c"));
    
    // detach b and add it back as a sub group of of "sub"
    flds->DetachGroup("b");
    flds->CreateGroup("sub")->AttachGroup(gb);
    
    EXPECT_TRUE(flds->HasGroup("a"));
    EXPECT_TRUE(flds->HasGroup("sub"));
    EXPECT_TRUE(flds->HasGroup("c"));

    EXPECT_EQ(flds->GetGroup("sub")->GetGroup("b"),gb);

    delete ds;
}






