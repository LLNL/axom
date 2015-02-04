//
// testfunc.cpp
//

#include <stdio.h>

#include "DatastoreInterface.hpp"
#include "DataFunctionCompiled.hpp"

void *work1(void)
{
    printf("In work1\n");
    return NULL;
}

int main(int argc, char *argv[])
{
  // Create DataStore
  DataStore::DataGroup* const myDS1 = DataStore::CreateDataStore("myDS1");

  DataStore::DataFunctionCompiled work1func(work1);
  DataStore::DataObject *dswork1 = myDS1->CreateDataObject("work1")->SetDataPointer(&work1func);

  work1func.Call();


  auto *dswork1_0 =  myDS1->GetData<DataStore::DataFunctionCompiled*>("work1");
  dswork1_0->Call();

  return 0;
}
