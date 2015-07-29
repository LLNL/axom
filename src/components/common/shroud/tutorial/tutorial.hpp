//
// tutorial.hpp
//

#include <string>

namespace tutorial
{

void Function1(void);
double Function2(double arg1, int arg2);
bool Function3(bool arg);
const std::string& Function4(const std::string& arg1, const std::string& arg2);
const std::string& Function4b(const std::string& arg1, const std::string& arg2);

class Class1
{
public:
    void Method1() {};
};

}  // namespace
