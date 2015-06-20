! C code that will be inserted into Python module via shroud splicer blocks

// splicer begin include
#include "sidre/sidre.hpp"
// splicer end include


// splicer begin class.DataStore.C_object
DataStore *ds;
// splicer end class.DataStore.C_object

// splicer begin class.DataStore.type.init
    PyErr_SetString(PyExc_NotImplementedError, "init");
    return -1;
// splicer end class.DataStore.type.init
