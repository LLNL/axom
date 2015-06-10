*******
2 Names
*******

This section contains guidelines for naming files, types (i.e., classes, 
structs, etc.), functions, data members, variables, etc. The main goal is 
to use a distinct and consistent naming convention for each item category 
so that the role of each entity in the code is obvious from the form of 
its name.
 
===========
2.1 General
===========

2.1.1 Every name **must** be meaningful. That is, its meaning **must** be clear to other code developers and users, not just the author of the name.

      A substantial benefit of good name selection is that it can greatly 
      reduce the amount of source code documentation required for others to 
      understand it. For example, when the name of a function clearly indicates 
      what it does and the meaning and purpose of each argument is clear from 
      its name, then code comments may be unnecessary. In fact, documenting 
      "the obvious" is often more distracting than helpful.

2.1.2 Each name **must** be consistent with other similar names in the code. 

      For example, if getter/setter methods follow the convention "getFoo" 
      and "setFoo" respectively, then adding a new setter method called 
      "putBar" is clearly inconsistent.

2.1.3 Tersely abbreviated or cryptic or names **should** be avoided. However, common acronyms and jargon that are well understood by team members **may** be used.


==============
2.2 File names
==============

In this section and throughout the guidelines, we refer to "header" and 
"source" files. Header files are those that are included in other files 
(either header or source). Source files are not included in other files and 
form translation units when a program is compiled. Also, header files often 
do not contain implementations, only function declarations and type definitions.Common exceptions include class template header files and headers that contain 
inline functions. Such files must contain implementations.

2.2.1 C++ header and source file extensions **must** be: \*.hpp and \*.cpp, respectively.

2.2.2 C-only header and source file extensions **must** be: \*.h and \*.c, respectively.

2.2.3 The name of each file **must** clearly indicate its contents.

      In partciular, a header file that defines a major type (e.g., a struct or 
      class) or contains function signatures associated with a major type 
      **must** include the type name in the file name. For example, the header 
      file for the class "MyClass" should be named "MyClass.hpp". 

      Similarly, a source file that implements functionality associated with 
      a major type, such as a class implementation file, **must** include 
      the type name in the file name. For example, the source file for the 
      class "MyClass" should be named "MyClass.cpp".

2.2.4 Header and source files that provide a limited, well-defined set of code functionality, but are not associated with a major data type, **must** be named so that the contents of the file are clear from the name.

      For example, header and source files containing unbound utility methods
      for dealing with file I/O **should** be named "FileUtils.hpp" and 
      "FileUtils.cpp", or similar.

2.2.5 File names that differ only in letter case **must** not be used.

      For example, having files with names "MyClass.hpp" and "myclass.hpp" 
      is not acceptable.


========================
2.3 Scope and type names
========================

2.3.1 All namespaces defined in the Toolkit **must** use all lowercase letters. 

      For example, most C++ entities in the CS Toolkit are inluded in the 
      namespace "asctoolkit"; for example::

         namespace asctoolkit {
              // . . .
         }

2.3.2 All type names (i.e., classes, structs, typedefs, enums, etc.) **must** be nouns and **must** be in mixed case with each word starting with an upper case
letter and all other letters in lower case.

      For example, these are acceptable type names::

         Line, BankAccount, MyNiftyClass

      These are not acceptable type names::

         line, LINE, Bankaccount, myNiftyClass

2.3.3 Separating characters, such as underscores, **should** not be used between words in a type name.

      For example, these are not acceptable type names::

         Bank_Account, My-Nifty-Class

2.3.4 Suffixes that may be used by compilers for name mangling, or which are used in the C or C++ standard libraries, such as "\_t", **must** not be used in type names.


==================
2.4 Function names
==================

2.4.1 Function names **must** use "camelCase " style.

      Specifically, the first word must be in lower case letters. If 
      multiple words are used, each word after the first must start with 
      an upper case letter and have all other letters in lower case. 
      Underscores must not be used in camelCase, but numbers may be 
      used. 

      For example, these are acceptable camelCase style function names::

         compAverage(), getName(), copy2()

2.4.2 All functions names (i.e., class members, unbound methods, etc.) **must** begin with a verb and clearly indicate what they do.

2.4.3 Complementary verbs (e.g., "get/set", "add/remove", "create/destroy") **must** be used for routines that perform complementary operations.

      Such symmetry prevents confusion and makes interfaces easier to use.

2.4.4 Verbs (e.g., "is", "has", "can") **must** be used for any function with a boolean return type (including an int value indicating true/false).

      For example, use::

         isInitialized(), hasLicense(), canEvaluate()


==================================
2.5 Data member and variable names
==================================

2.5.1 Variable names (e.g., class members, struct fields, function arguments, local variables) **must** use either "camelCase" style or "pot_hole" style. Within a Toolkit component, one consistent style **must** be used.

      **camelCase style:** The first word must be in lower case letters. 
      If multiple words are used, each word after the first must start with 
      an upper case letter and have all other letters in lower case. 
      Underscores must not be used in camelCase, but numbers may be used. 

      For example, these are acceptable camelCase style variable names::

         myAverage, personName, pressure2

      **pot_hole style:** All letters must be in lower case. If multiple 
      words are used, they must be separated by a single underscore. Numbers 
      may be used in pothole style names. 

      For example, these are acceptable pothole style variable names::

         my_average, person_name, pressure_2

2.5.2 Class and struct data member names **must** use one of the two prefixes: "m\_" and "s\_". 

      The prefix "m\_" indicates a regular data member and the prefix "s\_" 
      indicates a static member.

      This convention makes it obvious which variable names in the code refer 
      to class members/struct fields and which are local variables. For example,      the following are acceptable names for class data members using 
      camelCase style::

         m_myAverage, m_personName

      and pothole style::

         m_my_average, m_person_name

2.5.3 Verbs, such as "is", "has", "can", etc., **must** be used for each boolean variables, such as type bool or an integer that indicate true/false values.

      For example, use::

         m_is_initialized, has_license

      not::

         m_initialized, license

2.5.4 Local variables, such as loop indices, **should** be named so they are easy to search for using a text editor.

      For example, a loop index named "ivar" is easier to search for than 
      one named simply "i".

2.5.5 Each variable name **should** give an indication of its type.

      For example,::

         Topic* my_topic;

      is clearer than::

         Topic* my_value;


====================================
2.6 Macros and enumeration constants
====================================

2.6.1 Preprocessor macro constants **must** be named using all uppercase letters and underscores should be used between words.

      For example,::

         MAX_ITERATIONS, READ_MODE

2.6.2 The name of each enumeration value **should** start with a capital letter and use an underscore between words when multiple words are used.

       For example,::

          enum Orange
          {
             Navel,
             Valencia,
             Num_Orange_Types
          };

