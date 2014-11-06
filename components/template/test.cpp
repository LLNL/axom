
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "template.h"

TEST_CASE( "Objects can be constructed") {

    SampleClass sc = SampleClass{}; 

    SECTION(" Public functions behave as expected") { 
        REQUIRE( sc.returns1() == 1 ); 
        REQUIRE( sc.returns2() == 2); 
    }

    
}
