.. ##
.. ## Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory.
.. ##
.. ## LLNL-CODE-741217
.. ##
.. ## All rights reserved.
.. ##
.. ## This file is part of Axom.
.. ##
.. ## For details about use and distribution, please read axom/LICENSE.
.. ##

.. _portsec-label: 

===================================================
11 Portability, Compilation, and Dependencies
===================================================

C++ is a huge language with many advanced and powerful features. To avoid
over-indulgence and obfuscation, we would like to avoid C++ feature bloat.
By constraining or even banning the use of certain language features and
libraries, we hope to keep our code simple and portable. We also hope to 
avoid errors and problems that may occur when language features are not 
completely understood or not used consistently. This section lists such 
restrictions and explains why use of certain features is constrained or 
restricted.


--------------------------------------------------------------------
Portability
--------------------------------------------------------------------

Nothing beyond C++11
^^^^^^^^^^^^^^^^^^^^

11.1 C++ language features beyond standard C++11 **must not** be used unless
reviewed by the team and verified that the features are supported by all 
compilers we need to support.

      Changing this guideline requires full consensus of all team members.


Use C++11, but don't depend on it
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

11.2 Whenever C++11 features are used, an alternative implementation **must** 
be provided that conforms to the 2003 C++ standard.

      Applications that use Axom will rely on non-C++11 compilers 
      for our current generation of computing platforms, and possibly beyond, 
      so we must be able to compile and run our code with those compilers.
      Applications that use Axom will expect the code to compile
      and run with full functionality on all platforms they use. 

11.3 All C++11 usage **must** be guarded using the macro constant 
"AXOM_USE_CXX11" so that it can be compiled out of the code when necessary.

   For example, suppose you have a class that you want to support *move*
   semantics when available (i.e., when using a C++11-compliant compiler)
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
   #if defined(AXOM_USE_CXX11)
      /// Move constructor.
      MyClass(MyClass&& other);

      /// Move-assignment operator.
      MyClass& operator=(MyClass&& rhs);
   #endif

      // other methods...

   private:
      // data members...
   };


No non-standard language constructs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

11.4 Special non-standard language constructs, such as GNU extensions, 
**must not** be used if they hinder portability.


.. note:: Any deviation from these C++ usage requirements must be 
          agreed on by all members of the team and vetted with our
          main application users.


--------------------------------------------------------------------
Compilation
--------------------------------------------------------------------

Avoid conditional compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

11.5 Excessive use of the preprocessor for conditional compilation at a 
fine granularity (e.g., selectively including or removing individual source 
lines) **should** be avoided. 

      While it may seem convenient, this practice typically produces confusing 
      and error-prone code. Often, it is better to refactor the code into 
      separate routines or large code blocks subject to conditional compilation
      where it is more clear. 

      Code reviews by team members will dictate what is/is not acceptable.


The compiler is your friend
^^^^^^^^^^^^^^^^^^^^^^^^^^^

11.6 Developers **should** rely on compile-time and link-time errors to 
check for code correctness and invariants. 

      Errors that occur at run-time and which depend on specific control flow 
      and program state are inconvenient for users and can be difficult to 
      detect and fix.

.. important::  Add specific guidance on how this should be done...


--------------------------------------------------------------------
Minimize dependencies on third-party libraries (TPLs)
--------------------------------------------------------------------

11.7 While it is generally desirable to avoid recreating functionality that
others have already implemented, we **should** limit third-party library
dependencies for Axom to make it easier for users. We are a library,
and everything we necessarily depend on will become a dependency for our
user.  

      **Before introducing any significant TPL dependency on Axom
      (e.g., Boost), it must be agreed on by the development team and vetted
      with our main users.**

11.8 Unless absolutely necessary, any TPL we depend on **must not** be 
exposed through any public interface in Axom.
