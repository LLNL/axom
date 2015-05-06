#include "gtest/gtest.h"

#include "datastore/sidre.hpp"

using asctoolkit::sidre::DataBuffer;
using asctoolkit::sidre::DataGroup;
using asctoolkit::sidre::DataStore;
using asctoolkit::sidre::DataView;

using namespace conduit;

//------------------------------------------------------------------------------

TEST(datastore_view,create_views)
{
    DataStore *ds   = new DataStore();
    DataGroup *root = ds->getRoot();

    DataView *dv_0 = root->createViewAndBuffer("field0");
    DataView *dv_1 = root->createViewAndBuffer("field1");
    
    
    DataBuffer *db_0 = dv_0->getBuffer();
    DataBuffer *db_1 = dv_1->getBuffer();
        
    EXPECT_EQ(db_0->getUID(),0u);
    EXPECT_EQ(db_1->getUID(),1u);
    delete ds;
}

TEST(datastore_view,uint32_buffer_from_view)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    
    DataView *dv = root->createViewAndBuffer("u0");

    dv->allocate(DataType::uint32(10));
    uint32 *data_ptr = dv->getNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
        data_ptr[i] = i*i;

    dv->getNode().print_detailed();

    EXPECT_EQ(dv->getNode().schema().total_bytes(),
              dv->getDescriptor().total_bytes());
    delete ds;
    
}


TEST(datastore_view,uint32_array_multi_view)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataBuffer *dbuff = ds->createBuffer();

    dbuff->declare(DataType::uint32(10));
    dbuff->allocate();
    uint32 *data_ptr = dbuff->getNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
        data_ptr[i] = i;

    dbuff->getNode().print_detailed();

    EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
              dbuff->getDescriptor().total_bytes());


    DataView *dv_e = root->createView("even",dbuff);
    DataView *dv_o = root->createView("odd",dbuff);
    
    dv_e->apply(DataType::uint32(5,0,8));
    
    dv_o->apply(DataType::uint32(5,4,8));

    dv_e->getNode().print_detailed();
    dv_o->getNode().print_detailed();

    uint32_array dv_e_ptr = dv_e->getNode().as_uint32_array();
    uint32_array dv_o_ptr = dv_o->getNode().as_uint32_array();
    for(int i=0;i<5;i++)
    {
        std::cout << "idx:" <<  i  
                  << " e:" << dv_e_ptr[i] 
                  << " o:" << dv_o_ptr[i] 
                  << " em:" << dv_e_ptr[i]  % 2
                  << " om:" << dv_o_ptr[i]  % 2
                  << std::endl;

        EXPECT_EQ(dv_e_ptr[i] % 2,0u);
        EXPECT_EQ(dv_o_ptr[i] % 2,1u);
    }
    ds->print();
    delete ds;
    
}


TEST(datastore_view,init_uint32_array_multi_view)
{
    DataStore *ds = new DataStore();
    DataGroup *root = ds->getRoot();
    DataBuffer *dbuff = ds->createBuffer();

    dbuff->allocate(DataType::uint32(10));
    uint32 *data_ptr = dbuff->getNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
        data_ptr[i] = i;

    dbuff->getNode().print_detailed();

    EXPECT_EQ(dbuff->getNode().schema().total_bytes(),
              dbuff->getDescriptor().total_bytes());


    DataView *dv_e = root->createView("even",dbuff);
    DataView *dv_o = root->createView("odd",dbuff);
    
    // uint32(num_elems, offset, stride)
    dv_e->apply(DataType::uint32(5,0,8));


    // uint32(num_elems, offset, stride)    
    dv_o->apply(DataType::uint32(5,4,8));


    dv_e->getNode().print_detailed();
    dv_o->getNode().print_detailed();

    uint32_array dv_e_ptr = dv_e->getNode().as_uint32_array();
    uint32_array dv_o_ptr = dv_o->getNode().as_uint32_array();
    for(int i=0;i<5;i++)
    {
        std::cout << "idx:" <<  i  
                  << " e:" << dv_e_ptr[i] 
                  << " o:" << dv_o_ptr[i] 
                  << " em:" << dv_e_ptr[i]  % 2
                  << " om:" << dv_o_ptr[i]  % 2
                  << std::endl;

        EXPECT_EQ(dv_e_ptr[i] % 2,0u);
        EXPECT_EQ(dv_o_ptr[i] % 2,1u);
    }
    ds->print();
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
    DataGroup *root = ds->getRoot();
    
    // create a group to hold the "old" or data we want to copy
    DataGroup *r_old = root->createGroup("r_old");
    // create a view to hold the base buffer
    DataView  *base_old = r_old->createViewAndBuffer("base_data");

    // alloc our buffer
    // we will create 4 sub views of this array
    base_old->allocate(DataType::uint32(40));
    uint32 *data_ptr = base_old->getNode().as_uint32_ptr();
    
    
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
    DataBuffer *buff_old = base_old->getBuffer();
    buff_old->getNode().print();
    DataView *r0_old = r_old->createView("r0",buff_old);
    DataView *r1_old = r_old->createView("r1",buff_old);
    DataView *r2_old = r_old->createView("r2",buff_old);
    DataView *r3_old = r_old->createView("r3",buff_old);
    
    // each view is offset by 10 * the # of bytes in a uint32
    // uint32(num_elems, offset)
    index_t offset =0;
    r0_old->apply(DataType::uint32(10,offset));
    
    offset += sizeof(uint32) * 10;
    r1_old->apply(DataType::uint32(10,offset));
    
    offset += sizeof(uint32) * 10;
    r2_old->apply(DataType::uint32(10,offset));
    
    offset += sizeof(uint32) * 10;
    r3_old->apply(DataType::uint32(10,offset));

    /// check that our views actually point to the expected data
    //
    uint32 *r0_ptr = r0_old->getNode().as_uint32_ptr();
    for(int i=0;i<10;i++)
    { 
        EXPECT_EQ(r0_ptr[i],1u);
        // check pointer relation
        EXPECT_EQ(&r0_ptr[i],&data_ptr[i]);
    }
    
    uint32 *r3_ptr = r3_old->getNode().as_uint32_ptr();
    for(int i=0;i<10;i++)
    { 
        EXPECT_EQ(r3_ptr[i],4u);
        // check pointer relation
        EXPECT_EQ(&r3_ptr[i],&data_ptr[i+30]);
    }

    // create a group to hold the "old" or data we want to copy into
    DataGroup *r_new = root->createGroup("r_new");
    // create a view to hold the base buffer
    DataView  *base_new = r_new->createViewAndBuffer("base_data");

    // alloc our buffer
    // create a buffer to hold larger subarrays
    base_new->allocate(DataType::uint32(4 * 12));

    DataBuffer *buff_new = base_new->getBuffer();
    buff_new->getNode().print();

    // create the 4 sub views of this array
    DataView *r0_new = r_new->createView("r0",buff_new);
    DataView *r1_new = r_new->createView("r1",buff_new);
    DataView *r2_new = r_new->createView("r2",buff_new);
    DataView *r3_new = r_new->createView("r3",buff_new);
    
    // apply views to r0,r1,r2,r3
    // each view is offset by 12 * the # of bytes in a uint32

    // uint32(num_elems, offset)
    offset =0;
    r0_new->apply(DataType::uint32(12,offset));
    
    offset += sizeof(uint32) * 12;
    r1_new->apply(DataType::uint32(12,offset));
    
    offset += sizeof(uint32) * 12;
    r2_new->apply(DataType::uint32(12,offset));
    
    offset += sizeof(uint32) * 12;
    r3_new->apply(DataType::uint32(12,offset));

    /// update r2 as an example first
    buff_new->getNode().print();
    r2_new->getNode().print();
    
    /// copy the subset of value
    r2_new->getNode().update(r2_old->getNode());
    r2_new->getNode().print();
    buff_new->getNode().print();
    
    
    /// check pointer values
    uint32 *r2_new_ptr = r2_new->getNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
    { 
        EXPECT_EQ(r2_new_ptr[i],3u);
    }

    for(int i=10;i<12;i++)
    { 
        EXPECT_EQ(r2_new_ptr[i],0u); // assumes zero-ed alloc
    }


    /// update the other views
    r0_new->getNode().update(r0_old->getNode());
    r1_new->getNode().update(r1_old->getNode());
    r3_new->getNode().update(r3_old->getNode());
    
    buff_new->getNode().print();

    
    ds->print();
    delete ds;
    
}

TEST(datastore_view,simple_opaque)
{
    // create our main data store
    DataStore *ds = new DataStore();
    // get access to our root data Group
    DataGroup *root = ds->getRoot();
    int *src_data = new int[1];
    
    src_data[0] = 42;

    void *src_ptr = (void*)src_data;
    
    DataView *opq_view = root->createOpaqueView("my_opaque",src_ptr);
    
    // we shouldn't have any buffers
    EXPECT_EQ(ds->getNumberOfBuffers(),0u);
    
    EXPECT_TRUE(opq_view->isOpaque());
    
    void *opq_ptr = opq_view->getOpaque();
    
    int *out_data = (int*)opq_ptr;
    EXPECT_EQ(opq_ptr,src_ptr);
    EXPECT_EQ(out_data[0],42);
    
    ds->print();
    delete ds;
    delete [] src_data;
}




