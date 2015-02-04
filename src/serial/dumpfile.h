//
// dumpfile.h - Dump DataStore objects to a file.
//

#ifndef _DS_DUMPFILE_H
#define _DS_DUMPFILE_H

#include "pdb.h"
#include "DatastoreInterface.hpp"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct s_DF_DumpFile DF_DumpFile;

struct s_DF_DumpFile {
    PDBfile *fp;
    // error: size of array ‘type_map’ has non-integral type ‘DataStore::rtTypes::TypeID’
    //    char *type_map[DataStore::rtTypes::TypeID::undefined];
    char *type_map[9];
};


DF_DumpFile *DF_open(const char *name, const char *mode);
void DF_close(DF_DumpFile *fp);
void DF_write_group(DF_DumpFile *fp, DataStore::DataGroup *grp);


#ifdef __cplusplus
}
#endif


#endif
