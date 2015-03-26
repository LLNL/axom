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
    
    db0->SetDescriptor(DataType::uint64(10));
    db0->Allocate();
    db0->ApplyDescriptor();
    uint64 *db0_ptr = db0->GetNode().as_uint64_ptr();
    
    db1->SetDescriptor(DataType::float64(10));
    db1->Allocate();
    db1->ApplyDescriptor();
    float64 *db1_ptr = db1->GetNode().as_float64_ptr();
    
    for(int i=0;i<10;i++)
    {
        db0_ptr[i] = i;
        db1_ptr[i] = i*i;
    }
    
    DataView *f0 = flds->CreateView("u0",db0);
    DataView *f1 = flds->CreateView("f1",db1);
    
    
    flds->GetView("u0")->GetNode().print_detailed();
    flds->GetView("f1")->GetNode().print_detailed();
    ds->Print();
    delete ds;
}

