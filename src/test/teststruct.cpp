//
// teststruct.cpp - test structure types
//


#include <stdio.h>

#include "testds.hpp"

int main(int argc, char *argv[])
{
  DataStore::DataGroup* const myDS1 = DataStore::CreateDataStore("myDS1");

  test_add_struct1(myDS1);


}
