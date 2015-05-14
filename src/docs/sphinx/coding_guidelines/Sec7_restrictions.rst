************************************************
7 Restrictions on Language Usage and Libraries
************************************************

C++ is a huge language with many advanced and powerful features. To avoid
over-indulgence and obfuscation, we would like to avoid C++ feature bloat.
By constraining, or even banning, the use of certain language features and
libraries  we hope to keep code simple, portable, and avoid errors and 
problems that may occur when language features and usage are not completely 
understood or used consistently. This section lists such restrictions and 
explains why use of certain features is constrained or restricted.

==============
7.1 C++11
==============

C++11 may be used in limited ways as described in this section. Any other usage
must be carefully reviewed and approved by all team members.

7.1.1 All C++11 usage **must** be guarded using the macro constant "USE_CXX11" so that it can be compiled out of the code when necessary. 

      Applications that use the CS Toolkit will rely on non-C++11 compilers for
      our current generation of computing platforms, and possibly beyond, so we
      must be able to compile and run our code with those compilers.

7.1.2 Whenever C++11 features are used, an alternative implementation **must** be provided that conforms to the 2003 C++ standard.

      Applications that use the CS Toolkit will expect the code able to compile
      and run with full functionality on all platforms they use. 

**WE NEED TO WORK ON THIS**

==============
7.2 C++14
==============

No C++14 features are allowed at this time due to substantially incomplete 
compiler support on the platforms we care most about.

==============
7.2 Boost
==============

**WE NEED TO WORK ON THIS**
