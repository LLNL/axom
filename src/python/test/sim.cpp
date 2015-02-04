//
// sim.cpp - A sample simulation with embedded Python
//

#include <stdio.h>
#include "Python.h"
#include "numpy/arrayobject.h"

PyMODINIT_FUNC PyInit_datastore(void);
PyMODINIT_FUNC PyInit_dstest(void);

void start_python(void *root, const char *filename)
{
    int ierr;
    PyObject *m;
    PyObject *sim_m;
    PyObject *main_m, *main_d;

    Py_VerboseFlag = 0;
    //Py_SetProgramName("sim");
    //Py_SetPythonHome(Xstr(PYTHONHOME));
    ierr  = PyImport_AppendInittab("datastore", PyInit_datastore);
    ierr  = PyImport_AppendInittab("dstest", PyInit_dstest);
    Py_Initialize();

    main_m = PyImport_AddModule("__main__");
    main_d = PyModule_GetDict(main_m);

    /* Python: import datastore as ds */
    sim_m = PyImport_ImportModule("datastore");
    if (sim_m == nullptr) {
    }
    ierr = PyDict_SetItemString(main_d, "ds", sim_m);

    /* Python: import dstest */
    sim_m = PyImport_ImportModule("dstest");
    if (sim_m == nullptr) {
    }
    ierr = PyDict_SetItemString(main_d, "dstest", sim_m);

#if 0
    /* XXXX - Do not set excepthook yet */
    bas_d = PyModule_GetDict(bas_m);
    excepthook = PyDict_GetItemString(bas_d, "excepthook");
    ierr = PySys_SetObject("excepthook", excepthook);
#endif
    
    /* Python: import numpy as np */
    m = PyImport_ImportModule("numpy");
    if (m == nullptr) {
#if 0
	char *text = PDB_convert_exception();
	if (text != NULL) {
	    printf("%s\n", text);
	}
#endif
	printf("Unable to import numpy\n");
	exit(1);
    }
    ierr = PyDict_SetItemString(main_d, "np", m);

    if (_import_array() == -1) {
        printf("Unable to import numpy");
	exit(1);
    }

#if 0
    /* Add main module to code */
    if (root != NULL) {
	// XXX check for legal name (not "rdb" or "np")
	PyObject *obj = PDB_node_to_instance((TGT_node *) root);
	const char *name = TGT_node_get_name((TGT_node *) root);
	PyDict_SetItemString(main_d, name, obj);
    }
#endif

    if (filename != nullptr) {
	FILE *fp = fopen(filename, "r");
	if (fp != nullptr) {
	    PyRun_SimpleFile(fp, filename);
	}
    }

    PyRun_InteractiveLoop(stdin, "?");
    Py_Finalize();

    return;
}

//
//  usage:  exe [ file ]
//
int main(int argc, char *argv[])
{
    const char *filename = NULL;

    for (int i=1; i < argc; i++) {
	filename = argv[i];
    }

    start_python(nullptr, filename);
}
