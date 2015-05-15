// wrapDataGroup.h
// For C users and C++ implementation

#ifndef WRAPDATAGROUP_H
#define WRAPDATAGROUP_H

#include "wrapDataView.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef EXAMPLE_WRAPPER_IMPL
typedef void DS_datagroup;
#else
struct s_DS_datagroup;
typedef struct s_DS_datagroup DS_datagroup;
#endif

DS_dataview * DS_datagroup_create_view_and_buffer(DS_datagroup * self, const char * name);

DS_datagroup * DS_datagroup_create_group(DS_datagroup * self, const char * name);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATAGROUP_H
