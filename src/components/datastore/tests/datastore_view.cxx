#include "gtest/gtest.h"

#include "datastore/datastore.h"

using namespace DataStoreNS;
using namespace conduit;

//------------------------------------------------------------------------------

TEST(datastore_view,create_views)
{
    DataStore *ds   = new DataStore();
    DataGroup *root = ds->GetRoot();

    //
    // note we should prob have conv ways to create a view from the ds
    // 
    
    
    DataView *dv_0 = new DataView("field0",root,ds);
    DataView *dv_1 = new DataView("field1",root,ds);


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
    
    DataView *dv = root->CreateView("u0");

    dv->SetDescriptor(DataType::Arrays::uint32(10));
    dv->Allocate();
    dv->ApplyDescriptor();
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

    dbuff->SetDescriptor(DataType::Arrays::uint32(10));
    dbuff->Allocate();
    dbuff->ApplyDescriptor();
    uint32 *data_ptr = dbuff->GetNode().as_uint32_ptr();
    
    for(int i=0;i<10;i++)
        data_ptr[i] = i;

    dbuff->GetNode().print_detailed();

    EXPECT_EQ(dbuff->GetNode().schema().total_bytes(),
              dbuff->GetDescriptor().total_bytes());


    DataView *dv_e = new DataView("even",root,dbuff);
    DataView *dv_o = new DataView("odd",root,dbuff);
    
    dv_e->SetDescriptor(DataType::Arrays::uint32(5,0,8));
    dv_e->ApplyDescriptor();
    
    dv_o->SetDescriptor(DataType::Arrays::uint32(5,4,8));
    dv_o->ApplyDescriptor();

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

    delete ds;
    
}

