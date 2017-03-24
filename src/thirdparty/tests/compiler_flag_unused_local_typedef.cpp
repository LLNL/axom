#include <iostream>

/**
 * \file
 * This file tests the ability to disable warnings about unused local typedefs
 * on all supported compilers using the AXOM_DISABLE_UNUSED_LOCAL_TYPEDEF build variable.
 *
 */


int main()
{
    typedef int IntT;

    std::cout << "I have defined type IntT, but am not using it." << std::endl;
    return 0;
}


