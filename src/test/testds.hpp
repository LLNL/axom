//
// testds.hpp - Routines to build up datastore values for tests
//

#include "DatastoreInterface.hpp"
#if 1
#include "DataUserType.hpp"
#include "DataMember.hpp"
#endif

typedef struct s_struct1 struct1;

struct s_struct1 {
    int i1;
    int i2[10];
    int *i3;
};

void test_add_struct1(DataStore::DataGroup* grp);

