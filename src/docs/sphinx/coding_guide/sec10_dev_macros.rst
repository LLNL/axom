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

=======================================================
10 Code Development Macros
=======================================================

This section provides guidelines for consistent use of macros defined in the 
"common" Toolkit component and in our build system that we use in day-to-day 
code development.

.. important:: Fill in this section with guidelines for using the macros


------------------------------------
Error handling
------------------------------------

Runtime checks for incorrect or questionable usage and generation of 
informative messages that help users and developers use our code correctly.

"Sanity checks" should be performed on values of function arguments
(e.g., range checking, null pointer checking, etc.) upon entry to a function.

      This is an excellent way to provide run-time debugging capabilities in
      code. We have macros for this to make syntax consistent. When triggered,
      they can emit a failed boolean expression and descriptive message that
      help to understand the violation. They are active or not based on the
      compilation mode, either debug (active) or optimized (inactive). For
      example::

         void doSomething(int in_val, Foo* in_foo)
         {
            ATK_ASSERT_MSG( in_val >= 0, "in_val must be positive or zero" );
            ATK_ASSERT( in_foo != ATK_NULL_PTR );

            // ...function implementation...
         }


------------------------------------
Unused variables
------------------------------------

ATK_NOT_USED, etc.


------------------------------------
Disabling compiler generated methods
------------------------------------

DISABLE_COPY_AND_ASSIGNMENT


------------------------------------
Comditionally compiled code
------------------------------------





