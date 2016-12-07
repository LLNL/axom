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
