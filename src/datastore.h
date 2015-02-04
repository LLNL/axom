/*
 * datastore.h - C API for datastore.
 */

#ifndef __DS_DATASTORE_H
#define __DS_DATASTORE_H

#ifdef __cplusplus
extern "C" {
#endif

// In C++ DS_object is a void object to allow static_cast to work.
// In C it is an opaque object in an effort to help type safety.
// C may never dereference the pointer since is points to a C++ object.
#ifdef __cplusplus
typedef void DS_object;
#else
struct s_DS_object;
typedef struct s_DS_object DS_object;
#endif

int DS_init_library(void);
void DS_fin_library(void);


/* DATASTORE CREATION and EXTRACTION */

DS_object *DS_create_datastore(const char *name);

DS_object *DS_lookup_datastore(const char *name);

/* DATAGROUP CREATION and EXTRACTION */

DS_object *DS_create_datagroup(DS_object *dg, const char *name);

const char *DS_get_name(DS_object *obj);

/* DATA ALLOCATION/INSERTION */

/* create DS managed data object */

DS_object *DS_create_object(DS_object *ds, const char *name);

DS_object *DS_get_object(DS_object *ds, const char *name);  // DS_lookup_object

void *DS_allocate(DS_object *ds);

#if 0
DS_insert_array_into_datastore(const char *path, void *ptr);
DS_insert_array_into_group(DS_object *group, void *ptr);


DS_set_bool_attribute(DS_object *object, const char *name, bool value);
#endif

#ifdef __cplusplus
}
#endif

#endif
