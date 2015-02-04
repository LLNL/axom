//
// testdump.cpp
//

#include <stdio.h>

#include "DatastoreInterface.hpp"
#include "dumpfile.h"

int main(int argc, char *argv[])
{
  // Create DataStore
  DataStore::DataGroup* const problem = DataStore::CreateDataStore("myDS1");

  *(problem->CreateDataObject("numElems")
    ->SetType<int>()
    ->SetDataShape(1)
    ->Allocate()
    ->GetData<int*>()) = 10;
  *(problem->CreateDataObject("numFaces")
    ->SetType<int>()
    ->SetDataShape(1)
    ->Allocate()
    ->GetData<int*>()) = 20;

  DF_DumpFile *fp = DF_open("abc", "w");

  DF_write_group(fp, problem);

  DF_close(fp);

  return 0;
}
