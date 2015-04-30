#ifndef MESHAPI_UTILITIES_H_
#define MESHAPI_UTILITIES_H_

#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>

namespace {

    template<typename StringType>
    void killProcess(StringType const& msg)
    {
        std::cout << msg << std::endl;
        std::cout.flush();
        std::cerr.flush();
        abort();
    }

}


namespace asctoolkit{
namespace meshapi{

    typedef unsigned int        MeshIndexType;
    typedef unsigned long long  MeshSizeType;


    class NotImplementedException{};


// DBC assertions from Kull
#ifdef ATK_DEBUG
    #define DBC_ASSERTION(x, msg, kind) \
     if (!(x)) { \
      std::ostringstream s; \
      s << std::setprecision( 16 ) << std::setiosflags( std::ios_base::scientific ); \
      s << kind << ": (" << #x << ") is false: " << msg << '\n'; \
      s << "...at line " << __LINE__ << \
         " of file " << __FILE__ << "." << std::ends;\
      killProcess(s);\
     }

    #define REQUIRE2(x, msg) DBC_ASSERTION(x, msg, "Precondition violated")
    #define ASSERT2(x, msg) DBC_ASSERTION(x, msg, "Assertion violated")
    #define ENSURE2(x, msg) DBC_ASSERTION(x, msg,"Postcondition violated")
    #define INVARIANT2(x, msg) DBC_ASSERTION(x, msg, "Invariant violated")
#else
    #define ASSERT2(x, msg)
    #define REQUIRE2(x, msg)
    #define ENSURE2(x, msg)
    #define INVARIANT2(x, msg)
#endif

#ifdef ASSERT
#undef ASSERT
#endif
#define ASSERT(x) ASSERT2(x, "")
#define REQUIRE(x) REQUIRE2(x, "")
#define ENSURE(x) ENSURE2(x, "")
#define INVARIANT(x) INVARIANT2(x, "")



} // end namespace meshapi
} // end namespace asctoolkit

#endif //  MESHAPI_UTILITIES_H_
