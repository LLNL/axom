#include "Python.h"
//#include <pyTutorialmodule.hpp>
#include <stdio.h>

PyMODINIT_FUNC inittutorial(void);

int main(int argc, char** argv)  
{
    char filename[] = "test.py";
    FILE* fp;

    Py_Initialize();
    inittutorial();
    
    fp = fopen(filename, "r");
    PyRun_SimpleFile(fp, filename);
    fclose(fp);
    Py_Exit(0);  
    return 0;
}
