//
// ctest.cpp
//

#include "tutorial.hpp"


int main(int argc, char *argv[])
{
    tutorial::Function1();

    tutorial::Class1 *cptr = new tutorial::Class1();

    cptr->Method1();

    // Arguments
    // Integer and Real
    double rv2 = tutorial::Function2(1.5, 2);

    tutorial::Function6("name");
    tutorial::Function6(1);

}
