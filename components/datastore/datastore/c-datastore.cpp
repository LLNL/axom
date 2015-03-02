/*
 * datastore.c - C API for datastore.
 */

#include "datastore.h"
#include "DatastoreInterface.hpp"
#include "DataGroup.hpp"

extern "C" {

  DS_object *DS_create_datastore(const char *name)
  {
    DataStoreNS::DataGroup* const myDS1 = DataStoreNS::CreateDataStore(name);
    
    return (DS_object *) myDS1;
  }

  DS_object *DS_create_datagroup(DS_object *dg, const char *name)
  {
    DataStoreNS::DataGroup* group = static_cast<DataStoreNS::DataGroup *>( (void *) dg)->CreateDataGroup(name);
    return (DS_object *) group;
  }

  const char *DS_get_name(DS_object *obj)
  {
    return static_cast<DataStoreNS::DataObject *>( (void *) obj)->Name().c_str();
  }


}
