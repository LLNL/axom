// wrapDataStore.h
// For C users and C++ implementation

#ifndef WRAPDATASTORE_H
#define WRAPDATASTORE_H

#include "wrapDataBuffer.h"
#include "wrapDataGroup.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef EXAMPLE_WRAPPER_IMPL
typedef void DS_datastore;
#else
struct s_DS_datastore;
typedef struct s_DS_datastore DS_datastore;
#endif

DS_datastore * DS_datastore_new();

void DS_datastore_delete(DS_datastore * self);

DS_databuffer * DS_datastore_create_buffer(DS_datastore * self);

DS_datagroup * DS_datastore_get_root(DS_datastore * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPDATASTORE_H
