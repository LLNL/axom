//
//  test1c.c
//

#include "datastore.h"

int main(int argc, char *argv[])
{
  // Create DataStore
  DS_object *myDS1 = DS_create_datastore("myDS1");

  // use top level datastore group to create group
  DS_object *group1 = DS_create_datagroup(myDS1, "group1");

  const char *name = DS_get_name(group1);

  return 0;
}
