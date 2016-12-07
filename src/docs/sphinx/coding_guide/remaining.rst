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


=====================================
7 Code Formatting
=====================================

--------------------------------------------------------------------
7.1 Conditional statements and loops
--------------------------------------------------------------------

7.1.1 Curly braces **should** be used in all conditionals, loops, etc. 
even when the content inside the braces is a "one-liner". 

       This helps prevent coding errors and misinterpretation of intent. 
       For example, this::

          if (done) { ... }

       is preferable to this::

          if (done) ...

7.1.2 One-liners **may** be used for "if" conditionals with 
"else/else if"  clauses when the resulting code is clear. 

       For example, either of the following styles **may** be used::

          if (done) {
             id = 3;
          } else {
             id = 0;
          }

       or::

          if (done) { id = 3; } else { id = 0; }

7.1.3 Complex "if/else if" conditionals with many "else if" clauses 
**should** be avoided.

      Such statements can always be refactored using local boolean variables 
      or "switch" statements. Doing so often makes code easier to read and 
      understand and may improve performance.

7.1.4 An explicit test for zero/nonzero **must** be used in a conditional 
unless the tested quantity is a boolean or pointer type. 

      For example, a conditional based on an integer value should use::

         if (num_lines != 0) {

      not::

         if (num_lines) {


--------------------------------------------------------------------
7.2 White Space and Code Alignment
--------------------------------------------------------------------

Most conventions for indentation, spacing and code alignment 
preferred by the team are enforced by using the `uncrustify` tool. 
It can be run from the top-level CS Toolkit directory...

.. note :: Show how to run uncrustify on the code and where the format
           options are defined.

Not all preferred formatting conventions are supported by uncrustify.
The following guidelines provide additional recommendations to make
code easier to read and understand.


Use White Space To Make Code Easier to Read
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.2.1 Blank lines and indentation **should** be used throughout code to 
enhance readability. 

      Examples of helpful white space include:

         * Between operands and operators in arithmetic expressions.
         * After reserved words, such as "while", "for", "if", "switch", etc. 
           and before the parenthesis or curly brace that follows.
         * After commas separating arguments in functions.
         * After semi-colons in for-loop expressions.
         * Before and after curly braces in almost all cases.

7.2.2 White space **must not** appear between a function name and the opening 
parenthesis to the argument list. In particular, if a function call is broken 
across source lines, the break **must not** come between the function name and 
the opening parenthesis.

7.2.3 Tabs **must not** be used for indentation since this can be problematic 
for developers with different text editor settings.


Align Code Vertically to Show Scope
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.2.4 When function arguments (in either a declaration or implementation)
appear on multiple lines, the arguments **should** be aligned vertically 
for readability.

7.2.5 All statements within a function body **should** be indented within the surrounding curly braces.

7.2.6 All source lines in the same scope **should** be aligned vertically.
Continuation of previous lines **may** be indented if it make the code easier
to read.


Break Lines Where It Makes Sense
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.2.7 When a line is broken at a comma or semi-colon, it **must** be broken 
after the comma or semi-colon, not before. 

7.2.8 When a source line is broken at an arithmetic operator 
(i.e., , -, etc.), it **should** be broken after the operator, not before. 


Use Parentheses For Clarity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

7.2.9 Parentheses **should** be used in non-trivial mathematical and logical 
expressions to clearly indicate structure and intended order of operations. 
Do not assume everyone who reads the code knows all the rules for operator 
precedence.


.. _portsec-label: 

===================================================
8 Portability, Compilation, and Dependencies
===================================================

C++ is a huge language with many advanced and powerful features. To avoid
over-indulgence and obfuscation, we would like to avoid C++ feature bloat.
By constraining, or even banning, the use of certain language features and
libraries we hope to keep code simple, portable, and avoid errors and
problems that may occur when language features are not completely
understood or used consistently. This section lists such restrictions and
explains why use of certain features is constrained or restricted.


--------------------------------------------------------------------
8.1 Portability
--------------------------------------------------------------------

8.1.1 C++ language features beyond standard C++11 **must not** be used unless
reviewed by the team and verified that the features are supported by all 
compilers we need to support.

      Changing this guideline requires full con census of all team members.

8.1.2 Whenever C++11 features are used, an alternative implementation **must** 
be provided that conforms to the 2003 C++ standard.

      Applications that use the CS Toolkit will rely on non-C++11 compilers 
      for our current generation of computing platforms, and possibly beyond, 
      so we must be able to compile and run our code with those compilers.
      Applications that use the CS Toolkit will expect the code able to compile
      and run with full functionality on all platforms they use. 

8.1.3 All C++11 usage **must** be guarded using the macro constant "USE_CXX11" 
so that it can be compiled out of the code when necessary.

   For example, suppose you have a class that you want to support *move*
   semantics when available (i.e., when using a C++11-compilant compiler)
   and fall back to copy semantics otherwise:

.. code-block:: cpp

   class MyClass
   {
   public:

      /// Default constructor.
      MyClass();

      /// Destructor.
      ~MyClass();
      /// Copy constructor.
      MyClass(const MyClass& other);
      /// Copy-assignment operator.
      MyClass& operator=(const MyClass& rhs);
   #if defined(USE_CXX11)
      /// Move constructor.
      MyClass(MyClass&& other);

      /// Move-assignment operator.
      MyClass& operator=(MyClass&& rhs);
   #endif

      // other methods...

   private:
      // data members...
   };

8.1.4 Special non-standard language constructs, such as GNU extensions, 
**must not** be used if they hinder portability.


.. note :: Any deviation from these C++ usage requirements must be 
           agreed on by all members of the team and vetted with our
           main application users.


--------------------------------------------------------------------
8.2 Compilation
--------------------------------------------------------------------

8.2.1 Excessive use of the preprocessor for conditional compilation at a 
fine granularity (e.g., selectively including or removing individual source 
lines) **should** be avoided. 

      While it may seem convenient, this practice typically produces confusing 
      and error-prone code. Often, it is better to refactor the code into 
      separate routines or large code blocks subject to conditional compilation
      where it is more clear. 

      Code reviews by team members will dictate what is/is not acceptable.

8.2.2 Developers **should** rely on compile-time and link-time errors to 
check for code correctness and invariants. 

      Errors that occur at run-time and which depend on specific control flow 
      and program state are inconvenient for users and can be difficult to 
      detect and fix.

.. note ::      Add some specific guidance here on how this should be done...


--------------------------------------------------------------------
8.3 Third-party Library (TPL) Dependencies
--------------------------------------------------------------------

8.3.1 While it is generally desirable to avoid recreating functionality that
others have already implemented, we **should** limit third-party library
dependencies for the CS Toolkit to make it easier for users. We are a library,
and everything we necessarily depend on will become a dependency for our
user.  

      **Before introducing any significant TPL dependency on the Toolkit
      (e.g., Boost), it must be agreed on by the development team and vetted
      with our main users.**

8.3.2 Unless absolutely necessary, any TPL we depend on **must not** be 
exposed through any public interface in the CS Toolkit.


.. _codingrefs-label:

======================================
References and Useful Resources
======================================

Most of the guidelines here were gathered from the following list sources. 
The list contains a variety of useful resources for programming in C++
beyond what is presented in these guidelines.

#. *The Chromium Projects: C++ Dos and Don'ts*. https://www.chromium.org/developers/coding-style/cpp-dos-and-donts

#. Dewhurst, S., *C++ Gotchas: Avoiding Common Problems in Coding and Design*, Addison-Wesley, 2003.

#. Dewhurst S., *C++ Common Knowledge: Essential Intermediate Programming*, Addison-Wesley, 2005.

#. *Doxygen manual*, http://www.stack.nl/~dimitri/doxygen/manual/index.html

#. *Google C++ Style Guide*, https://google-styleguide.googlecode.com/svn/trunk/cppguide.html

#. *ISO/IEC 14882:2011 C++ Programming Language Standard*.

#. Josuttis, N., *The C++ Standard Library: A Tutorial and Reference, Second Edition*, Addison-Wesley, 2012.

#. Meyers, S., *More Effective C++: 35 New Ways to Improve Your Programs and Designs*, Addison-Wesley, 1996.

#. Meyers, S., *Effective STL: 50 Specific Ways to Improve Your Use of the Standard Template Library*, Addison-Wesley, 2001.

#. Meyers, S., *Effective C++: 55 Specific Ways to Improve Your Programs and Designs (3rd Edition)*, Addison-Wesley, 2005.

#. Meyers, S., *Effective Modern C++: 42 Specific Ways to Improve Your Use of C++11 and C++14*, O'Reilly.

#. Programming Research Ltd., *High-integrity C++ Coding Standard, Version 4.0*, 2013.

#. Sutter, H. and A. Alexandrescu, *C++ Coding Standards: 101 Rules, Guidelines, and Best Practices*, Addison-Wesley, 2005.
