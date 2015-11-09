/*
 *  ftutorial.c - C routines called from Fortran.
 */

#include "chasm/chasm_Fortran_binding.h"


void *global_allocatable;  // silence 'unused parameter' warnings

// called directly from Fortran XXX - must mangle name properly
void all_test1_(void *array)
{
  global_allocatable = array;
  CHASM_printArrayDesc(array, 1);
}
