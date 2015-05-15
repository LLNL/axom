// wrapDataStore.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapDataStore.h"
#include "sidre/DataStore.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

DS_datastore * DS_datastore_new()
{
DataStore *selfobj = new DataStore();
// splicer begin
return (DS_datastore *) selfobj;
// splicer end
}

void DS_datastore_delete(DS_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin
delete selfobj;
// splicer end
}

DS_databuffer * DS_datastore_create_buffer(DS_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin
DataBuffer * rv = selfobj->createBuffer();
return rv;
// splicer end
}

DS_datagroup * DS_datastore_get_root(DS_datastore * self)
{
DataStore *selfobj = static_cast<DataStore *>(self);
// splicer begin
DataGroup * rv = selfobj->getRoot();
return rv;
// splicer end
}

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"
