// wrapDataGroup.cpp
#define EXAMPLE_WRAPPER_IMPL
#include "wrapDataGroup.h"
#include "sidre/DataGroup.hpp"

extern "C" {
namespace asctoolkit {
namespace sidre {

DS_dataview * DS_datagroup_create_view_and_buffer(DS_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataView * rv = selfobj->createViewAndBuffer(name);
return rv;
// splicer end
}

DS_datagroup * DS_datagroup_create_group(DS_datagroup * self, const char * name)
{
DataGroup *selfobj = static_cast<DataGroup *>(self);
// splicer begin
DataGroup * rv = selfobj->createGroup(name);
return rv;
// splicer end
}

}  // namespace asctoolkit
}  // namespace sidre
}  // extern "C"
