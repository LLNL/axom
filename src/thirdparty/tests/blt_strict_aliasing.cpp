#include<iostream>

struct Foo {int i;  Foo* ob;};
struct Bar {int i;  Foo* ob;};

int main()
{
    Foo foo = { 1, NULL};
    ((Bar*)(&foo))->i++;      // violates strict aliasing

    std::cout << " foo.i: " << foo.i << std::endl;
    return 0;
}


