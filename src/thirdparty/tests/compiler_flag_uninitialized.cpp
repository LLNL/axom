#include<iostream>
#include<cstdlib>

/**
 * \file
 * This file tests the ability to disable warnings about uninitialized variables
 * on all supported compilers using the AXOM_DISABLE_UNINITIALIZED_WARNINGS build variable.
 */


int main()
{
    int* result;        // Note: variable not allocated or initialized

    if( rand()%2 == 0 )
      *result = 5;

    return 0;
}
