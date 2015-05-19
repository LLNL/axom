************************************************
7 Restrictions on Language Usage and Libraries
************************************************

C++ is a huge language with many advanced and powerful features. To avoid
over-indulgence and obfuscation, we would like to avoid C++ feature bloat.
By constraining, or even banning, the use of certain language features and
libraries we hope to keep code simple, portable, and avoid errors and 
problems that may occur when language features are not completely 
understood or used consistently.  This section lists such restrictions and 
explains why use of certain features is constrained or restricted.


=======================
7.1 C++11 and beyond
=======================

Applications that use the CS Toolkit will rely on non-C++11 compilers for 
our current generation of computing platforms, and possibly beyond, so we
must be able to compile and run our code with those compilers.

C++11 may be used in the CS Toolkit in limited ways as described in this 
section. Any other usage must be carefully reviewed and approved by all 
team members.

7.1.1 All C++11 usage **must** be guarded using the macro constant "USE_CXX11" so that it can be compiled out of the code when necessary. 

.. code-block:: cpp

   #if defined(USE_CXX11)
   #include <unordered_map>
   #else
   #include <boost/unordered_map>
   #endif

   // ...

   #if defined(USE_CXX11)
      typedef std::unordered_map<std::string, common::IDType> MapType;
   #else
      typedef boost::unordered_map<std::string, common::IDType> MapType;
   #endif

7.1.2 Whenever C++11 features are used, an alternative implementation **must** be provided that conforms to the 2003 C++ standard.

      Applications that use the CS Toolkit will expect the code able to compile
      and run with full functionality on all platforms they use. 

7.1.3 C++14 features **must not** be used due to substantially incomplete compiler support on the platforms we care most about.


**WE NEED TO WORK ON THIS SECTION**


==============
7.2 Boost
==============

The Boost C++ libraries are generally high quality and provide many powerful
and useful capabilities not found in the core C++ language. Indeed, some Boost
libraries eventually make their way into the C++ standard.
 
Some LLNL codes have used Boost successfully for many years. However, 
version inconsistencies (e.g., changes from one version of Boost to the next or 
two codes using different incompatible versions that need to be compiled
and linked into the same executable) and compiler portability have presented 
problems in the past. To avoid increasing the maintenace burden for 
applications that use the CS Toolkit, we restrict Boost usage in the CS 
Toolkit as described in this section. Any other usage must be carefully 
reviewed and approved by all team members.

7.2.1 All CS Toolkit components **must** use the same version of Boost that is maintained for the Toolkit.

7.2.2 Boost libraries that require compilation **must not** be used. That is, only those libraries that provide header files **may** be used.

7.2.3 Boost usage **must not** be exposed through any public interface in the CS Toolkit. 
