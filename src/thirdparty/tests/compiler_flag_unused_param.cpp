#include<iostream>

/**
 * \file
 * This file tests the ability to disable warnings about unused parameters
 * on all supported compilers using the ATK_DISABLE_UNUSED_PARAMETER_WARNINGS build variable.
 */

void foo(int param)
{
    std::cout << "Hello " << std::endl;
}


int main()
{
    foo(2);
    return 0;
}
