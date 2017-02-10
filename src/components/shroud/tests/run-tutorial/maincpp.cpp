//
// Tests for tutorial.cpp
//

#include "tutorial.hpp"


int main(int argc, char *argv[])
{
    tutorial::Class1 * obj = new tutorial::Class1;

    obj->Method1();

    delete obj;
}
