! C code that will be inserted into Python module via shroud splicer blocks

// splicer beginX include
#include "sidre/sidre.hpp"
// splicer endX include


// splicerX begin class.DataStore.C_object
DataStore *ds;
// splicerX end class.DataStore.C_object

// splicerXbegin class.DataStore.type.init
    PyErr_SetString(PyExc_NotImplementedError, "init");
    return -1;
// splicerXend class.DataStore.type.init
