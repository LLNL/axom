.. ##
.. ## Copyright (c) 2016, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory.
.. ##
.. ## All rights reserved.
.. ##
.. ## This file cannot be distributed without permission and
.. ## further review from Lawrence Livermore National Laboratory.
.. ##

.. _testing-label:

****************************************
Software Testing 
****************************************

This section describes how to write and run tests and perform other 
code health tasks.

===================
Writing Unit Tests
===================

Describe basic goals and requirements for unit tests...

Describe how to write and add unit tests...

===========================
Writing Integration Tests
===========================

Describe basic goals and requirements for integration tests...

Describe how to write and add integration tests...


============================
Google Test Tips
============================

We are using the macro-based gtest unit testing framework for the toolkit.
Simple tests

For simple unit tests, a test function is defined by the TEST macro
TEST(<test_group>,<test_name>)
{
 // ... code here
}
Some useful macros

 

Within each test, macros are used to define the expected values and test assertions.

    Failure of expected values will eventually cause the test to fail, but will allow it to continue executing.
    Failure of assertions will also cause the test to stop running. 


// testing for boolean values
EXPECT_TRUE( <condition>); 
EXPECT_FALSE( <condition>); 
ASSERT_TRUE( <condition>); 
ASSERT_FALSE( <condition>); 
 
// testing for equality
EXPECT_EQ ( <expected_value>, <actual_value>);
...
 
// fuzzy comparisons for floating point numbers, with and without explicit tolerance
EXPECT_FLOAT_EQ(expected, actual); 
EXPECT_DOUBLE_EQ(expected, actual); 
EXPECT_NEAR(val1, val2, abs_error);
 
// when you expect the code to crash
ASSERT_DEATH( { <code_or_function_call>, <regexp_of_failure> }


Note for the ASSERT_DEATH:  to avoid a warning about thread safety, you can add:

 
::testing::FLAGS_gtest_death_test_style = "threadsafe";

For more information, see: https://code.google.com/p/googletest/wiki/AdvancedGuide#Death_Tests
Failure messages

The gtest TEST macros overload the ostream operator<<, so you can add custom messages to any test.

E.g.
EXPECT_TRUE( predicate(val)) << " This test was supposed to pass with val = " << val ;
Test fixtures

Test fixtures require their own class derived from::testing::Test, which has virtual functions SetUp() and TearDown().

Tests that use this call the macro TEST_F(){}

See https://code.google.com/p/googletest/wiki/V1_7_Primer#Test_Fixtures:_Using_the_Same_Data_Configuration_for_Multiple_Te
Stub for main function

gtest will generate a main function for you if it is not explicitly defined.

However, if you need to write your own, e.g. for initialization code or pre-/post-processing, the following stub (which initializes the SLIC logging library) might be useful:
Stub for a gtest main function with SLIC initialization
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;
 
int main(int argc, char * argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);
  UnitTestLogger logger;  // create & initialize test logger,
                          // finalized when exiting main scope
  result = RUN_ALL_TESTS();
  return result;
}

 
Color in gtest console output

If you would like to have the gtest output its messages in color (e.g. green for pass and red for fail), you can set the GTEST_COLOR environment variable to 'yes'.  E.g. in bash:
export GTEST_COLOR=yes

 
Links

https://code.google.com/p/googletest/wiki/V1_7_Primer

https://code.google.com/p/googletest/wiki/AdvancedGuide


================================================
Filename and CMake Target Conventions for Tests
================================================

Note that following these conventions should clearly group all targets for each 
component together. Having consistent suffices at the end of each name should ma
ke it clear what the target is for.

    Test file names and targets
        We should also add a target for each component to run all its tests; e.g
., 'make sidre_tests'
        Test file name format:
            <component name>_<test name>_<optional language specifier>_test
        Examples:
            sidre_buffer_test.cpp (C++ code test, language specifier optional)
            sidre_buffer_C_test.cpp (C language test)
            sidre_buffer_F_test.f (Fortran language test).
        When named like this and added to the CMakeLists.txt file, the test targ
ets would appear somewhere in the list of make targets as follows when typing 'm
ake help':
            sidre_buffer_test
            sidre_buffer_C_test
            sidre_buffer_F_test
    Example file names and targets:
        We should also add a target for each component to run all its examples; 
e.g., 'make sidre_examples'
        Example file name format:
            <component name>_<example name>_<optional language specifier>_exampl
e
        Examples:
            sidre_shocktube_example.cpp (C++ example)
            sidre_shocktube_F_example.f (Fortran example)
        How they appear in 'make help':
            sidre_shocktube_example
            sidre_shocktube_F_example
    Documentation make targets
        Component user guide docs target name:
            <component name>_user_docs
        Component source code docs (i.e., doxygen) target name:
            <component name>_doxygen_docs
        For example:
            sidre_user_docs
            sidre_doxygen_docs
        Other documentation examples:
            coding_guide_docs
            dev_guide_docs
            quickstart_guide_docs
        When typing 'make help', each of these examples appears the same as its target name.


===============
Running Tests
===============

Describe how to run tests from the command line using make targets....

See :ref:`bamboo-label` for information about setting up and running tests 
using the *Bamboo* continuous integration (CI) tool.

