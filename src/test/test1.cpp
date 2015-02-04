//
// test1.cpp
//

#include <stdio.h>

#include "DatastoreInterface.hpp"

int main(int argc, char *argv[])
{
  // Create DataStore
  DataStore::DataGroup* const myDS1 = DataStore::CreateDataStore("myDS1");

  // get reference to a data store object.
  DataStore::DataGroup* const myDS1_2 = DataStore::GetDataStore("myDS1");

  //--------------------------------------------------

  // use top level datastore group to create group
  myDS1->CreateDataGroup("group1");

  // get reference to the group
  DataStore::DataGroup* const group1 = myDS1->GetDataGroup("group1");

  // use top level datastore group to create group, and assign to local reference
  DataStore::DataGroup* const group2 = myDS1->CreateDataGroup("group2");

  // use group to create group, and assign to local reference
  DataStore::DataGroup* const group1a = group1->CreateDataGroup("group1a");


  // iterate over groups

  for( auto obj : myDS1->GetDataObjects() ) {
    fprintf(stdout, "%s\n", obj->Name().c_str());
  }


  return 0;
}
