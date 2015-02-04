//
// testds.cpp - Routines to build up datastore values for tests
//

#include "testds.hpp"

void test_add_struct1(DataStore::DataGroup* grp)
{
    auto *type = new DataStore::DataUserType("struct1");

    type->AddMember("i1");
    type->AddMember("i2");

#if 0
    type->AddMember("i1", "int");
    type->AddMember("i2", "int", "[10]");
    type->AddMember("i3", "int*", "[i1]");


    inttype = LookupType("int")
    type->AddMember("i1", inttype);


    for (auto mem : type->GetMembers() ) {
    }
#endif


    return;
}

