//
// tutorial.hpp
//

#include <string>

namespace tutorial
{

void Function1(void);
double Function2(double arg1, int arg2);
bool Function3(bool arg);
const std::string& Function4a(const std::string& arg1, const std::string& arg2);
const std::string& Function4b(const std::string& arg1, const std::string& arg2);

double Function5(double arg1 = 3.13, int arg2 = 5);

void Function6(const std::string &name);
void Function6(int indx);

template<typename ArgType>
void Function7(ArgType arg)
{
    return;
}

template<typename RetType>
RetType Function8()
{
    return 0;
}

void Function9(double arg);

class Class1
{
public:
    void Method1() {};
};

}  // namespace
