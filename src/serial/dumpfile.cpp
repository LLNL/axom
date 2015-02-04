//
// dumpfile.c - Dump DataStore objects to a file.
//

#include <assert.h>
#include "dumpfile.h"

#ifdef __cplusplus
extern "C" {
#endif

DF_DumpFile *DF_open(const char *name, const char *mode)
{
    //    DF_DumpFile *fp = (DF_DumpFile *) malloc(sizeof(DF_DumpFile));
    DF_DumpFile *fp = new(DF_DumpFile);
    assert(fp != NULL);
    fp->fp = PD_open(name, const_cast<char *>(mode));
    assert(fp->fp != NULL);

    if (sizeof(int) == 4) {
	fp->type_map[static_cast<int>(DataStore::rtTypes::TypeID::int32_id)] = SC_strsave("int");
    }
    // XXX other types and sizes

    return fp;
}

void DF_close(DF_DumpFile *fp)
{
    assert(fp->fp != NULL);
    int ierr = PD_close(fp->fp);
    assert(ierr == TRUE);
    //    free(fp);
    delete fp;
    return;
}

void DF_write_group(DF_DumpFile *fp, DataStore::DataGroup *grp)
{
    for( auto obj : grp->GetDataObjects() ) {
	const char *name = obj->Name().c_str();
	char *type = fp->type_map[static_cast<int>(obj->GetType())];
	void *data = obj->GetData<void*>();
	long ind[21];

	ind[0] = 0;  // min
	ind[1] = obj->length() - 1; // max
	ind[2] = 1;  // stride

	int ok = PD_write_alt(fp->fp, name, type, data,
			      1, ind);
    }
}

#ifdef __cplusplus
}
#endif
