/*
 * datastore.c - C API for datastore.
 */

#include "datastore.h"
#include "DatastoreInterface.hpp"
#include "DataGroup.hpp"

extern "C" {

  DS_object *DS_create_datastore(const char *name)
  {
    sidre::DataGroup* const myDS1 = sidre::CreateDataStore(name);
    
    return (DS_object *) myDS1;
  }

  DS_object *DS_create_datagroup(DS_object *dg, const char *name)
  {
    sidre::DataGroup* group = static_cast<sidre::DataGroup *>( (void *) dg)->CreateDataGroup(name);
    return (DS_object *) group;
  }

  const char *DS_get_name(DS_object *obj)
  {
    return static_cast<sidre::DataObject *>( (void *) obj)->Name().c_str();
  }


}
