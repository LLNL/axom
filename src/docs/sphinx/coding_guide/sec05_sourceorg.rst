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

.. _sourceguide-label:

=====================================
5 Source File Organization
=====================================

The goal is to make it easy to find essential information in a file easily 
and quickly. Consistently-applied conventions for file organization
can significantly improve user understanding and developer productivity. 

---------------------------------------------------------
Each source file should have an associated header file
---------------------------------------------------------

5.1 Each source file **should** have an associated header file with a matching
name, such as "Foo.hpp" for the source file "Foo.cpp".

.. note:: **Exceptions:** Test files may not require headers.


---------------------------------------------------------
Avoid extraneous header file inclusions
---------------------------------------------------------

5.2 Unnecessary header files **should not** be included in source files 
(i.e., headers not needed to compile the file standalone).

      Such header file inclusions introduce spurious file dependencies, which
      may increases compilation time unnecessarily.


---------------------------------------------------------
Function order in source and header files should match
---------------------------------------------------------

5.3 The order of functions implemented in a source file **should** match the 
order in which they appear in the associated header file.

      This makes the methods easier to locate and compare with documentation
      in the header file.


.. _sourcelayout-label:

---------------------------------------------------------
Source file layout details
---------------------------------------------------------

Content **must** be organized consistently in all source files.
This section summarizes the recommended source file layout using numbers
and text to illustrate the basic structure. Details about individual items
are contained in the guidelines after the summary.

.. code-block:: cpp

   // (1) Doxygen file prologue

   // (2) CS Toolkit copyright and release statement

   // (3) Header file inclusions (only those that are NECESSARY)
   #include "..."

   // (4a) Toolkit project namespace declaration
   namespace asctoolkit {

   // (5a) Internal namespace (if used); e.g.,
   namespace awesome {

   // (6) Initialization of static variables and data members, if any; e.g.,
   Foo* MyClass::s_shared_foo = 0;

   // (7) Implementation of static class member functions, if any

   // (8) Implementation of non-static class members and other methods

   // (5b) Internal namespace closing brace (if needed)
   } // awesome namespace closing brace

   // (4b) Project namespace closing brace
   } // asctoolkit namespace closing brace


5.4 **(Item 1)** Each source file **must** begin with a Doxygen file prologue.

      See :ref:`docsec-label` for details.

5.5 **(Item 2)** Each source file **must** contain a comment section that 
includes the CS Toolkit copyright and release statement.

      See :ref:`docsec-label` for details.

5.6 **(Item 3)** All necessary header file include statements **must** appear 
immediately after the copyright and release statement and before any 
implementation statements in the file.

.. note :: If a header is included in a header file, it **should not** be 
           included in the associated source file.

5.7 **(Items 4a, 4b, 5a, 5b)** All contents in a source file **must** follow 
the same namespace inclusion pattern as its corresponding header file 
(See :ref:`headerlayout-label`).

      Either the main project namespace (item 4a) or internal namespace 
      (item 5a) may be used, or both may be used. A closing brace ( "}" ) 
      is required to close each namespace declaration (items 4b and 5b).

5.8 **(Item 6)** Any static variables and class data members that are 
defined in a header file **must** be initialized in the associated source 
file before any method implementations.

5.9 **(Items 7, 8)** Static method implementations **must** appear before 
non-static method implementations.
