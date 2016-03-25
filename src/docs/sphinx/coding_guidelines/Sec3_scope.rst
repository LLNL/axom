*******************************
3 Directories, Files, and Scope
*******************************

This section contains some basic directory and file organization guidelines. 
While each toolkit component may establish its own directory structure, 
following these guidelines will help make it easy to locate a specific file 
and, once the file is found, to locate essential information in it easily 
and quickly.


=====================================
3.1 File location and directory names
=====================================

It is common practice for C++ libraries to have header files and associated
implementation files located in the same directory. We follow this practice.

3.1.1 Each file **must** reside in a directory that corresponds to the code functionality supported by the contents of the file.

3.1.2 Each source directory **must** be named so that the collective purpose of the files it contains is clear. Each directory name **must** be in all lower case letters and should consist of a single word with no non-alphabetic characters.


============================================================
3.2 C-only header file (\*.h extension) content organization
============================================================

This section contains guidelines for C-only header files (with "\*.h" 
extension).  The file layout summary here uses numbers and text to illustrate 
the basic structure. Details about individual items follow.

.. code-block:: cpp

   /* (1) Copyright statement */

   /* (2) Doxygen file prologue */

   /* (3a) Header file include guard, e.g., */
   #ifndef MYCOOLSTUFF_H
   #define MYCOOLSTUFF_H

   /* (4) Header file inclusions (when NEEDED in lieu of forward declarations) */
   #include "..."

   /* (5) Forward declarations NEEDED in the header file */
   struct ...;

   /* (6) Type definitions (struct, enum, etc.) with Doxygen comments e.g., */

   /*!
    * \brief Brief ...summary comment text...
    *
    * ...detailed comment text...
    */
   typedef struct MyStructType {
      ...
   } MyStruct;

   /* (7a) Function prototypes declared with external linkage (if needed to support interaction with other packages) */
   #ifdef __cplusplus
   extern "C" {
   #endif

   /* (8) Prototypes for "extern" routines, e.g.,
   extern double MyCoolStuff_execFoo1(MyStuff_t* stuff);
   extern void MyCoolStuff_execFoo2(int* vec, int vec_len);

   /* (7b) External linkage closing brace (if needed) */
   #ifdef __cplusplus
   }
   #endif

   /* (3b) Header file include guard closing endif */
   #endif /* closing endif for header file include guard */


The numbers in parentheses in the following guidelines correspond to the 
numbered items in the comments in the preceding summary.

3.2.1 Each header file **must** begin with a comment section containing the LLNL copyright statement (item 1 in summary). 

      See Section 4 for details.
  
3.2.2 A Doxygen file prologue (item 2 in summary) **should** follow the copyright statement. 

      See Section 4 for details.

3.2.3 The contents of each header file **must** be guarded using a preprocessor directive that defines a unique "guard name" for the header file. 

      The guard must appear immediately after the file prologue and use the 
      #ifndef directive (item 3a in summary); this requires a closing #endif 
      statement at the end of the file (item 3b in summary). The preprocessor 
      constant must use the file name followed by "_H"; e.g., "MYCOOLSTUFF_H" 
      as above.

3.2.4 All necessary header file include statements (item 4 in summary) **must** appear immediately after definition of the header guard and before any forward
declarations, type definitions, or function prototypes.

3.2.5 Forward declarations (item 5 in summary) needed for type definitions and function prototypes **must** appear before all type definitions and prototypes.

3.2.6 All type definitions (item 6 in summary) **must** appear after the header file inclusions and before any subroutine prototypes. A proper type definition prologue **must** appear before each type definition; see Section 4 for details.

3.2.7 All function prototypes **must** be declared with external "C" linkage if and only if this is needed for compatibility with another software package (items 7a and 7b in summary).

      This is done by wrapping all prototypes within an "extern" statement, 
      which is guarded with a special macro constant::

         #ifdef __cplusplus
            extern "C" {
         #endif

            /* function prototypes... */

         #ifdef __cplusplus
         {
         #endif

3.2.8 Function prototypes **must** be declared using the "extern" keyword (item 8 in summary).


============================================================
3.3 C-only source file (\*.c extension) content organization
============================================================

Do we need this?


============================================================
3.4 C++ header file (\*.hpp extension) content organization
============================================================

This section contains guidelines for C++ header files (with "\*.hpp" extension).
The file layout summary here uses numbers and text to illustrate the basic 
structure. Details about individual items follow.

.. code-block:: cpp

   // (1) Copyright statement

   // (2) Doxygen file prologue

   // (3a) Header file include guard, e.g.,
   #ifndef MYCLASS_HPP
   #define MYCLASS_HPP

   // (4) Header file inclusions (when NEEDED in lieu of forward declarations)
   #include "..."

   // (5) Forward declarations NEEDED in header file (outside of asctoolkit
   namespace)
   class ...;

   // (6a) "asctoolkit" namespace declaration
   namespace asctoolkit {

   // (7a) Toolkit component namespace (if used); e.g.,
   namespace awesome {

   // (8) Forward declarations NEEDED in header file (in toolkit namespace(s)
   class ...;

   // (9) Type definitions (class, enum, etc.) with Doxygen comments e.g.,
   /*!
    * \brief Brief ...summary comment text...
    *
    * ...detailed comment text...
    */
   Class MyClass {
      ...
   } MyClass;

   // (7b) Toolkit component namespace closing brace (if needed)
   } // awesome namespace closing brace

   // (6b) "asctoolkit namespace closing brace
   } // asctoolkit namespace closing brace

   // (3b) Header file include guard closing endif */
   #endif // closing endif for header file include guard

The numbers in parentheses in the following guidelines correspond to the 
numbered items in the comments in the preceding summary.

3.4.1 Each header file **must** begin with a comment section containing the LLNL copyright statement (item 1 in summary). 

      See Section 4 for details.

3.4.2 A Doxygen file prologue (item 2 in summary) **should** follow the copyright statement. 

      See Section 4 for details.

3.4.3 The contents of each C++ header file **must** be guarded using a preprocessor directive that defines a unique "guard name" for the header file.

      The guard must appear immediately after the file prologue and use the 
      '#ifndef' directive (item 3a in summary); this requires a closing 
      '#endif' statement at the end of the file (item 3b in summary). The 
      preprocessor constant must use the file name followed by "_HPP"; e.g., 
      "MYCLASS_HPP" as above.

3.4.4 All necessary header file include statements (item 4 in summary) **must** appear immediately after definition of the header guard and before any forward
declarations, type definitions, etc.

3.4.5 Any necessary forward declarations (item 5 in summary) for types defined outside the toolkit namespace **must** appear before the toolkit namespace statement.

3.4.6 All types defined and methods defined in a C++ header file **must** be included in a namespace. 

      Either the main "asctoolkit" namespace (item 6a in summary) or a toolkit 
      component namespace (item 7a in summary) may be used, or both may be 
      used. A closing brace ( "}" ) is required to close each namespace 
      declaration (items 6b and 7b) before the closing '#endif' for the header 
      file include guard.

3.4.7 Forward declarations for types defined in the toolkit, and which are needed for the header file, **must** appear first in the "asctoolkit" or nested namespace before any other statements (item 8 in summary).

3.4.8 All class and other type definitions (item 9 in summary) **must** appear after the header file inclusions and forward declarations. A proper class prologue **must** appear before the class definition; see Section 4 for details.


============================================================
3.5 C++ source file (\*.cpp extension) content organization 
============================================================

This section contains guidelines for C++ source files (with "\*.cpp" extension).
The file layout summary here uses numbers and text to illustrate the basic 
structure. Details about individual items follow.

.. code-block:: cpp

   // (1) Copyright statement

   // (2) Header file inclusions (only those that are NECESSARY)
   #include "..."

   // (3a) "asctoolkit" namespace declaration
   namespace asctoolkit {

   // (4a) Toolkit component namespace (if used); e.g.,
   namespace awesome {

   // (5) Initialization of static class data members, if any; e.g.,
   Foo* MyClass::s_shared_foo = 0;

   // (6) Implementation of static class member functions, if any

   // (7) Implementation of non-static class members and other methods

   // (4b) Toolkit component namespace closing brace (if needed)
   } // awesome namespace closing brace

   // (3b) "asctoolkit namespace closing brace
   } // asctoolkit namespace closing brace

The numbers in parentheses in the following guidelines correspond to the 
numbered items in the comments in the preceding summary.

3.5.1 Each source file **must** begin with a comment section containing the LLNL copyright statement (item 1 in summary).

3.5.2 All necessary header file include statements (item 2 in summary) **must** appear immediately after the copyright statement and before any actual implementation statements in the file.

3.5.3 All contents in a C++ source file **must** follow the same namespace inclusion pattern as its corresponding header file (see item 3.4.6). 

      Either the main "asctoolkit" namespace (item 3a in summary) or a toolkit 
      component namespace (item 4a in summary) may be used, or both may be used.
      A closing brace ( "}" ) is required to close each namespace declaration 
      (items 3b and 4b) before the closing '#endif' for the header file include
      guard.

3.5.4 When used, static class data members **must** be initialized explicitly 
in the class source file before any member functions are defined (item 5 in summary).

3.5.6 All implementations of static class member functions (item 6 in summary), if any, **must** appear before implementations of non-static class member functions (item 7 in summary).


==================================
3.6 General header file guidelines
==================================

Good header file structure and conventions can make a huge positive impact on 
readability, and productivity of software developers. In earlier sections, we 
described basic header file organizational guidelines. In this section, we 
provide additional header file guidelines.

3.6.1 Each source file **must** have an associated header file with a matching name, such as "Foo.hpp" for the source file Foo.cpp".

      **Exceptions:** Unit test files and the file containing main do not 
      require headers.

3.6.2 Header files **may** contain multiple type definitions (e.g., structs, classes, enums, etc.). However, type definitions and function declarations in a header file **must** be related closely and/or support the primary type for which the file is named.

3.6.3 A header file **must** be self-contained and self-sufficient.

      In particular, a header file

      * Must have proper header file include guards (as illustrated in previous         sections) to prevent multiple inclusion. The macro symbol name for each 
        guard must be chosen to guarantee uniqueness within a compilation unit.
      * Must include all other headers and/or forward declarations it needs to 
        be compiled (i.e., each type used in the header file must be accounted 
        for). In addition, a file should not rely on symbols defined in another 
        header file that it includes; the other file should be included 
        explicitly.
      * Must contain the implementations of all generic templates and inline 
        methods defined in it. A compiler will require the full definitions of 
        these constructs to be seen in every source file that uses them. 

        **Exceptions:** Function templates or class template members whose 
        implementations are fully specialized with all template arguments must 
        be defined in an associated source file to avoid linker errors. Fully 
        specialized templates are not templates and so they are treated just 
        like any other function.

3.6.4 Header files **should** use forward declarations instead of header file inclusions when possible.

      This avoids having the compiler open more files than are needed, which 
      can speed up recompilation when header files change.

      **Exceptions:** 

      * Header files that define external APIs for Toolkit components **must** 
        include all header files for all types that appear in the API. This 
        makes use of the API much easier.
      * When using a function, such as an inline method or template, that is 
        implemented in a header file, the header file containing the 
        implementation must be included.
      * Similarly, when using C+ standard library types in a header file, it 
        **may** be preferable to include the associated standard headers in the 
        header file to make it easier to use. This avoids having explicit 
        inclusion of standard headers wherever the header file is used.

3.6.5 A forward type declaration **must** be used in a header file when an include statement would result in a circular dependency among header files or when the only the type name is needed and not the type definition.

3.6.6 Unnecessary header files or forward declarations (i.e., when a type definition or name is not needed) **should not** be included in header files.

      Such header file inclusions, in particular, introduce spurious file 
      dependencies, which unnecessarily increases the number of files that 
      are opened during code compilation.

3.6.7 Header file include statements **should** use the same ordering pattern for all files within a toolkit component. 

      This improves code readability, helps to avoid misunderstanding 
      dependencies, and insures successful compilation regardless of 
      dependencies in other files. A common header file inclusion ordering 
      scheme is:

      1. Related header (e.g., class header in class implementation file)
      2. C library headers
      3. C++ library headers
      4. Headers from other libraries
      5. Project headers

      Also, code is easier to understand when include files are ordered 
      alphabetically within each of these sections and a blank line is 
      inserted between sections. Also, adding comments that describe the 
      header file categories are sometimes useful.  For example,

.. code-block:: cpp

         // Related header
         #include "MyClass.hpp"

         // C standard library (including non-std unistd.h)
         #include <stdio.h>
         #include <unistd.h>

         // C++ standard library
         #include <hash_map>
         #include <vector>

         // "base" library headers
         #include "base/Port.hxx"

         // Headers from this project
         #include "MyOtherClass.hpp"

3.6.8 A "typedef" statement, defining a synonymous name for a type, **should** appear in the header file where the type is defined. In addition, a header file **should** only define a synonymous name for a type whose definition appears in that same header file.

      These practices help insure that all names associated with a given type 
      are available when the appropriate header file is used and eliminates 
      potentially inconsistent type names.

3.6.9 Routines **should** be ordered and grouped in a header file to enhance 
code readability and understanding.

      For example, all related methods should be grouped together.

3.6.10 The name of each function argument **must** be specified in a header file declaration. Also, names in function declarations and definitions **must** match.

       For example, this is not an acceptable function declaration::

          void doSomething(int, int, int);

3.6.11 Each function, type, and variable declaration in a header file **must** be documented according to the guidelines in Section 4.

       However, clear names that are self-explanatory are typically preferable 
       to reduce the need to write (and maintain!) documentation. For example,
       short, simple functions (e.g., inline functions) with related 
       functionality should be grouped together and described with a single 
       prologue if the resulting documentation is clearer and more concise.


==================================
3.7 General source file guidelines
==================================

3.7.1 Unnecessary header files **should not** be included in source files (i.e., not needed to compile the file).

      Such header file inclusions introduce spurious file dependencies, which 
      unnecessarily increases the number of files that are opened during code 
      compilation.

3.7.2 The order of routines implemented in a source file **should** match the order in which they appear in the associated header file.

      This makes the methods easier to locate and compare with documentation 
      in the header file.

3.7.3 Each function implementation in a source file **should** be documented according to the guidelines in Section 4.


==========
3.8 Scope
==========

3.8.1 All C++ code in the toolkit **must** be included in a namespace. 

      Either the main "asctoolkit" namespace or a toolkit component namespace 
      **may** be used, or both **may** be used with the component namespace 
      nested within the "asctoolkit" namespace.

3.8.2 When a toolkit component namespace is used, it **must** be unique within the toolkit.

      In particular, Toolkit components **must** not share a namespace.

3.8.3 The C++ using directive **must not** be used in any header file.  

      Using this directive in a header file leverages a bad decision to 
      circumvent the namespace across every file that directly or indirectly 
      includes that header file. Note that this guideline implies that each 
      type name appearing in a header file **must be fully-qualified** (i.e., 
      using the namespace identifier and scope operator) if it resides in a 
      different namespace than the contents of the file.

3.8.4 The C++ using directive **may** be used in source files to avoid the need to use a fully-qualified type name at each declaration. When used, using directives **must** appear after all "#include" directives in the file.

3.8.5 When only parts of a namespace are used in an implementation file, only those parts **should** be included with a using directive instead of the entire namespace contents.

      For example, if you only need the standard library vector container form 
      the "std" namespace, it is preferable to use::

         using std::vector;

      rather than::

         using namespace std;

3.8.6 Non-member functions that are meant to be used only in a single source file **should** be placed in the unnamed namespace to limit their scope to that file.

      This guarantees link-time name conflicts will not occur. For example::

         namespace {
            void myInternalFunction();
         }

3.8.7 Nested classes **should** be private unless they are part of the enclosing class interface.

      For example::

         class Outer
         {
            // ...
         private:
            class Inner
            {
               // ...
            };
         };

      When only the enclosing class uses a nested class, making it private 
      does not pollute the outer scope needlessly. Furthermore, nested classes 
      can be forward declared within the enclosing class definition and then 
      defined in the implementation file for the enclosing class. For example::

         class Outer
         {
            class Inner; // forward declaration

            // use name 'Inner' in Outer class definition
         };

         // In Outer.cxx implementation file...
         class Outer::Inner
         {
            // Inner class definition
         }

      This makes it clear that the nested class is only needed in the 
      implementation and does not clutter the class definition.

3.8.8 Local variables **should** be declared in the narrowest scope possible and as close to first use as possible.

      Minimizing variable scope makes source code easier to comprehend and 
      may also have performance benefits. For example, declaring a loop index 
      inside a for-loop statement such as::

         for (int ii = 0; ...) {

      is preferable to::

         int ii;
         ...
         for (ii = 0; ...) {

      **Exception:** When a local variable is an object, its constructor and 
      destructor may be invoked every time a scope (such as a loop) is entered 
      and exited, respectively. Thus, instead of this::

         for (int ii = 0; ii < 1000000; ++ii) {
            Foo f;
            f.doSomethingCool(ii);
         }

      it may be more efficient to do this::

         Foo f;
         for (int ii = 0; ii < 1000000; ++ii) {
            f.doSomethingCool(ii);
         }

3.8.9 Static or global variables of class type **must not** be used. 

      Due to indeterminate order of construction, their use may cause bugs 
      that are very hard to find. Static or global variables that are pointers 
      to class types **may** be used and must be initialized properly in a 
      single source file.

3.8.10 A reference to any item in the global namespace (which should be rare if needed at all) **should** use the scope operator ("::") to make this clear.

      For example::

         int local_val = ::global;
