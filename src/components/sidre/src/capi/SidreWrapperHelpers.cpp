/**
 *  \file SidreWrapperHelpers.hpp
 *
 *  \brief File used to contain helper functions for Fortran/C API wrappers.
 *         User code should not include this file.
 *
 */

#include "SidreWrapperHelpers.hpp"

namespace asctoolkit
{
namespace sidre
{

static char *global_char;
static int global_int;
static void *global_void;

void *register_allocatable(char *name, int lname,
			   void *array, int atk_type, int rank)
{
    global_char = name;
    global_int = lname;
    global_void = array;
    global_int = atk_type;
    global_int = rank;
    return NULL;
}

} /* end namespace sidre */
} /* end namespace asctoolkit */

