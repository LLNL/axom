#include "gtest/gtest.h"

#include "datastore/datastore.h"

using namespace DataStoreNS;
using namespace conduit;

//------------------------------------------------------------------------------

TEST(datastore_view,create_views)
{
    DataStore *ds   = new DataStore();
    DataGroup *root = ds->GetRoot();

    DataView *dv_0 = root->CreateViewAndBuffer("field0");
    DataView *dv_1 = root->CreateViewAndBuffer("field1");
    
    
    DataBuffer *db_0 = dv_0->GetBuffer();
    DataBuffer *db_1 = dv_1->GetBuffer();
        
    EXPECT_EQ(db_0->GetUID(),0);
    EXPECT_EQ(db_1->GetUID(),1);
    delete ds;
}

TEST(datastore_view,uint32_buffer_from_view)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->GetRoot();
    
    DataView *dv = root->CreateViewAndBuffer("u0");

    dv->Allocate(DataType::uint32(10));
    uint32 *data_ptr = dv->GetNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
        data_ptr[i] = i*i;

    dv->GetNode().print_detailed();

    EXPECT_EQ(dv->GetNode().schema().total_bytes(),
              dv->GetDescriptor().total_bytes());
    delete ds;
    
}


TEST(datastore_view,uint32_array_multi_view)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->GetRoot();
    DataBuffer *dbuff = ds->CreateBuffer();

    dbuff->Declare(DataType::uint32(10));
    dbuff->Allocate();
    uint32 *data_ptr = dbuff->GetNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
        data_ptr[i] = i;

    dbuff->GetNode().print_detailed();

    EXPECT_EQ(dbuff->GetNode().schema().total_bytes(),
              dbuff->GetDescriptor().total_bytes());


    DataView *dv_e = root->CreateView("even",dbuff);
    DataView *dv_o = root->CreateView("odd",dbuff);
    
    dv_e->Apply(DataType::uint32(5,0,8));
    
    dv_o->Apply(DataType::uint32(5,4,8));

    dv_e->GetNode().print_detailed();
    dv_o->GetNode().print_detailed();

    uint32_array dv_e_ptr = dv_e->GetNode().as_uint32_array();
    uint32_array dv_o_ptr = dv_o->GetNode().as_uint32_array();
    for(int i=0;i<5;i++)
    {
        std::cout << "idx:" <<  i  
                  << " e:" << dv_e_ptr[i] 
                  << " o:" << dv_o_ptr[i] 
                  << " em:" << dv_e_ptr[i]  % 2
                  << " om:" << dv_o_ptr[i]  % 2
                  << std::endl;

        EXPECT_EQ(dv_e_ptr[i] % 2,0);
        EXPECT_EQ(dv_o_ptr[i] % 2,1);
    }
    ds->Print();
    delete ds;
    
}


TEST(datastore_view,init_uint32_array_multi_view)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->GetRoot();
    DataBuffer *dbuff = ds->CreateBuffer();

    dbuff->Allocate(DataType::uint32(10));
    uint32 *data_ptr = dbuff->GetNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
        data_ptr[i] = i;

    dbuff->GetNode().print_detailed();

    EXPECT_EQ(dbuff->GetNode().schema().total_bytes(),
              dbuff->GetDescriptor().total_bytes());


    DataView *dv_e = root->CreateView("even",dbuff);
    DataView *dv_o = root->CreateView("odd",dbuff);
    
    // uint32(num_elems, offset, stride)
    dv_e->Apply(DataType::uint32(5,0,8));


    // uint32(num_elems, offset, stride)    
    dv_o->Apply(DataType::uint32(5,4,8));


    dv_e->GetNode().print_detailed();
    dv_o->GetNode().print_detailed();

    uint32_array dv_e_ptr = dv_e->GetNode().as_uint32_array();
    uint32_array dv_o_ptr = dv_o->GetNode().as_uint32_array();
    for(int i=0;i<5;i++)
    {
        std::cout << "idx:" <<  i  
                  << " e:" << dv_e_ptr[i] 
                  << " o:" << dv_o_ptr[i] 
                  << " em:" << dv_e_ptr[i]  % 2
                  << " om:" << dv_o_ptr[i]  % 2
                  << std::endl;

        EXPECT_EQ(dv_e_ptr[i] % 2,0);
        EXPECT_EQ(dv_o_ptr[i] % 2,1);
    }
    ds->Print();
    delete ds;
    
}


TEST(datastore_view,uint32_array_multi_view_resize)
{
    ///
    /// This example creates a 4 * 10 buffer of ints,
    /// and 4 views that point the 4 sections of 10 ints
    ///
    /// We then create a new buffer to support 4*12 ints
    /// and 4 views that point into them
    ///
    /// after this we use the old buffers to copy the values
    /// into the new views
    ///
    
    // create our main data store
    DataStore *ds = new DataStore();
    // get access to our root data Group
    DataGroup *root = ds->GetRoot();
    
    // create a group to hold the "old" or data we want to copy
    DataGroup *r_old = root->CreateGroup("r_old");
    // create a view to hold the base buffer
    DataView  *base_old = r_old->CreateViewAndBuffer("base_data");

    // alloc our buffer
    // we will create 4 sub views of this array
    base_old->Allocate(DataType::uint32(40));
    uint32 *data_ptr = base_old->GetNode().as_uint32_ptr();
    
    
    // init the buff with values that align with the
    // 4 subsections.
    for(int i=0;i<10;i++)
        data_ptr[i] = 1;
    for(int i=10;i<20;i++)
        data_ptr[i] = 2;
    for(int i=20;i<30;i++)
        data_ptr[i] = 3;
    for(int i=30;i<40;i++)
        data_ptr[i] = 4;


    /// setup our 4 views
    DataBuffer *buff_old = base_old->GetBuffer();
    buff_old->GetNode().print();
    DataView *r0_old = r_old->CreateView("r0",buff_old);
    DataView *r1_old = r_old->CreateView("r1",buff_old);
    DataView *r2_old = r_old->CreateView("r2",buff_old);
    DataView *r3_old = r_old->CreateView("r3",buff_old);
    
    // each view is offset by 10 * the # of bytes in a uint32
    // uint32(num_elems, offset)
    index_t offset =0;
    r0_old->Apply(DataType::uint32(10,offset));
    
    offset += sizeof(uint32) * 10;
    r1_old->Apply(DataType::uint32(10,offset));
    
    offset += sizeof(uint32) * 10;
    r2_old->Apply(DataType::uint32(10,offset));
    
    offset += sizeof(uint32) * 10;
    r3_old->Apply(DataType::uint32(10,offset));

    /// check that our views actually point to the expected data
    //
    uint32 *r0_ptr = r0_old->GetNode().as_uint32_ptr();
    for(int i=0;i<10;i++)
    { 
        EXPECT_EQ(r0_ptr[i],1);
        // check pointer relation
        EXPECT_EQ(&r0_ptr[i],&data_ptr[i]);
    }
    
    uint32 *r3_ptr = r3_old->GetNode().as_uint32_ptr();
    for(int i=0;i<10;i++)
    { 
        EXPECT_EQ(r3_ptr[i],4);
        // check pointer relation
        EXPECT_EQ(&r3_ptr[i],&data_ptr[i+30]);
    }

    // create a group to hold the "old" or data we want to copy into
    DataGroup *r_new = root->CreateGroup("r_new");
    // create a view to hold the base buffer
    DataView  *base_new = r_new->CreateViewAndBuffer("base_data");

    // alloc our buffer
    // create a buffer to hold larger subarrays
    base_new->Allocate(DataType::uint32(4 * 12));

    DataBuffer *buff_new = base_new->GetBuffer();
    buff_new->GetNode().print();

    // create the 4 sub views of this array
    DataView *r0_new = r_new->CreateView("r0",buff_new);
    DataView *r1_new = r_new->CreateView("r1",buff_new);
    DataView *r2_new = r_new->CreateView("r2",buff_new);
    DataView *r3_new = r_new->CreateView("r3",buff_new);
    
    // apply views to r0,r1,r2,r3
    // each view is offset by 12 * the # of bytes in a uint32

    // uint32(num_elems, offset)
    offset =0;
    r0_new->Apply(DataType::uint32(12,offset));
    
    offset += sizeof(uint32) * 12;
    r1_new->Apply(DataType::uint32(12,offset));
    
    offset += sizeof(uint32) * 12;
    r2_new->Apply(DataType::uint32(12,offset));
    
    offset += sizeof(uint32) * 12;
    r3_new->Apply(DataType::uint32(12,offset));

    /// update r2 as an example first
    buff_new->GetNode().print();
    r2_new->GetNode().print();
    
    /// copy the subset of value
    r2_new->GetNode().update(r2_old->GetNode());
    r2_new->GetNode().print();
    buff_new->GetNode().print();
    
    
    /// check pointer values
    uint32 *r2_new_ptr = r2_new->GetNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
    { 
        EXPECT_EQ(r2_new_ptr[i],3);
    }

    for(int i=10;i<12;i++)
    { 
        EXPECT_EQ(r2_new_ptr[i],0); // assumes zero-ed alloc
    }


    /// update the other views
    r0_new->GetNode().update(r0_old->GetNode());
    r1_new->GetNode().update(r1_old->GetNode());
    r3_new->GetNode().update(r3_old->GetNode());
    
    buff_new->GetNode().print();

    
    ds->Print();
    delete ds;
    
}

/// TODO: Add Tests for 
/// *opaque
/// CreateViewAndBuffer and DestroyViewAndBuffer
/// DestroyView, checking buffer state and the state of other views






