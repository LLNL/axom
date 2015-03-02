#include "Catch.hpp"

#include "template.hpp"

SCENARIO( "A Sample Class is pretty silly", "[template]") { 
    GIVEN ("A Stack allocated Sample Class Object") { 
        SampleClass sc = SampleClass{}; 
        WHEN("we call some methods"){
            THEN(" they return what they say they should") { 
                REQUIRE(sc.returns1() == 1); 
                REQUIRE(sc.returns2() == 2); 
            }
        }
    }
}
