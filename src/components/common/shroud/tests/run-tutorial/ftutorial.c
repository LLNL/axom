/*
 *  ftutorial.c - C routines called from Fortran.
 */

//#include "compilers/GNU.h"
// F90_SetCompilerCharacteristics()
int printArrayDesc_GNU(const void* desc, int rank);


void *global_allocatable;  // silence 'unused parameter' warnings

// called directly from Fortran XXX - must mangle name properly
void all_test1_(void *array)
{
  global_allocatable = array;
  printArrayDesc_GNU(array, 1);
}
