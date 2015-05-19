****************************
5 General Code Development
****************************

This section contains various development guidelines intended to improve code 
readability, correctness, portability, consistency, and robustness.


=====================================================
5.1 General design and implementation considerations
=====================================================

5.1.1 Simplicity, clarity, ease of modification and extension **should** always be a main goal when writing new code or changing existing code. 

5.1.2 All designs and implementations **should** be reviewed with other team members and refined based on input from others. 

      This is especially important for designs that are complex or potentially 
      unclear. What cannot be easily understood cannot be changed and 
      maintained with confidence.

5.1.3 Each entity (class, struct, variable, function, etc.) **should** embody one clear, well-defined concept. 

      The responsibilities of an entity may increase as it is used in new and 
      different ways. However, changes that divert it from its original intent 
      **should** be avoided. Also, large, monolithic entities that provide too 
      much functionality or which include too many concepts tend to increase 
      code coupling and complexity and introduce undesirable side effects. 
      Smaller, clearly constrained objects are easier to write, test, maintain,
      and use correctly. Also, small, simple objects tend to get used more often
      and reduce code redundancy. Designs and implementations that are overly 
      complexity should be evaluated by the team and modified appropriately.

5.1.4 Global, complex, or opaque data sharing **should** be avoided. Shared data increases coupling and contention between different parts of a code base, which makes maintenance and modification difficult.

5.1.5 When making substantial modifications or stylistic changes to existing code, an attempt **should** be made to make all other code, for example in a source file, consistent with the changes.


=====================
5.2 Code robustness 
=====================

5.2.1 The "const" qualifier **should** be used for variables and methods when appropriate to clearly indicate usage and to take advantage of compiler-based error-checking. 

      Constant declarations make code safer and less error-prone since they 
      enforce intent at compile time. They also simplify code understanding
      because a constant declaration clearly indicates the fact that the state
      of a variable or object will not change in the scope in which the 
      declaration appears.

5.2.2 Preprocessor macros **should not** be used when there is a better alternative, such as an inline function or a constant variable definition. 

      For example, this::

         const double PI = 3.1415926535897932384626433832;

      is preferable to this::

         #define PI (3.1415926535897932384626433832)

      Macros circumvent the ability of a compiler to enforce beneficial 
      language concepts such as scope and type safety. Macros are also 
      context-specific and can produce errors that cannot be understood 
      easily in a debugger. Macros **should be used only** when they are the 
      best choice for a particular situation.

5.2.3 An enumeration type **should** be used instead of macro definitions or "int" data for sets of related constant values. 

      In C++, enums are distinct types with a compile-time specified set of 
      values. Enumeration values cannot be implicitly cast to integers or 
      vice versa -- a "static_cast" operator must be used to make the 
      conversion explicit. Thus, enums provide type and value safety and 
      scoping benefits.

5.2.4 Hard-coded numerical constants and other "magic numbers" **must not** be used directly in code. When such values are needed, they **should** be declared as named constants to enhance code readability and consistency.

5.2.5 Floating point constants **should** always be written with a decimal point and have at least one digit before and after the decimal point for clarity. 

      For example, use "0.5" instead of ".5" and "1.0" instead of "1" or "1.". 


=================================
5.3 Compilation and portability
=================================

5.3.1 All C-only files **must** contain only standard C99 usage. Use of standard C11 features **must** be agreed upon by the project team and be guarded in the code using the "USE_C11" compiler generated macro constant. 

      Changing this guideline requires full concensus of all team members.

5.3.2 All C++ files **must** contain only standard C++03 usage. Use of standard C++11 or C++14 features **must** be agreed upon by the project team. If C++11 standard features are introduced, they **must** be guarded in the code using the "USE_CXX11" compiler generated macro constant. 

      Changing this guideline requires full concensus of all team members.

5.3.3 Special non-standard language constructs, such as GNU extensions, **must not** be used if they hinder portability.

5.3.4 Excessive use of the preprocessor for conditional compilation at a fine granularity (e.g., selectively including or removing individual source lines) **should** be avoided. 

      While it may seem convenient, this practice typically produces confusing 
      and error-prone code. Often, it is better to refactor the code into 
      separate routines or large code blocks subject to conditional compilation
      where it is obvious. The team **should** establish a policy policy for 
      how this is done.

5.3.5 Developers **should** rely on compile-time and link-time errors to check for code correctness and invariants. 

      Errors that occur at run-time and which depend on specific control flow 
      and program state are inconvenient for users and can be difficult to 
      detect and fix.

5.3.6 Before committing code to the source repository, developers **must** attempt to compile cleanly at the highest warning level with the main compiler(s) supported by the project. All warnings **must** be understood and eliminated if possible (not by reducing the warning level!). 

      Compiler warnings, while seemingly innocuous at times, often indicate 
      problems that do not appear until later or until specific run-time 
      conditions are encountered.


=======================
5.4 Memory management
=======================

5.4.1 Memory **should** be deallocated in the same scope in which it is allocated.

5.4.2 Memory **should** be deallocated as soon as it is no longer needed.

5.4.3 Pointers **should** be set to null explicitly when memory is deallocated. 

      For uniformity across the CS Toolkit and to facilitate C++11 and 
      non-C++11 usage, this should be done using the common macro 
      "ATK\_NULLPTR"; For example:: 

         double* data = new double[10];
         // ...
         delete [ ] data;
         data = ATK_NULLPTR;
  
5.4.4 Data managed exclusively within C++ code **must** be allocated and deallocated using the "new" and "delete" operators. 

      The operator "new" is type-safe, simpler to use, and less error-prone 
      than the "malloc" family of C functions.  C++ new/delete operators 
      **must not** be combined with C malloc/free functions.

5.4.5 Every C++ array deallocation statement **must** include "[ ]" (i.e., "delete[ ]") to avoid memory leaks. 

      The rule of thumb is: when "[ ]" appears in the allocation, then "[ ]" 
      **must** appear in the corresponding deallocation statement.  

5.4.6 Before committing code to the source repository, one **should** use memory-checking tools to verify there are no leaks and other memory misuse.

      When merging to the *develop* or *master* branches, compilation with a 
      variety fo compilers, testing, memory-checking, etc. will be done 
      automatically as part of the *pull request* approval process.  The pull
      request will not be approved until all of these tasks succeed.


===========================
5.5 Function declarations
===========================

5.5.1 Any class member function that does not change a data member of the associated class **must** be declared "const".

5.5.2 Function arguments **must** be ordered the same way for all routines in a project.

      Common conventions are either to put all input arguments first, then 
      outputs, or the other way around. Input and output and outputs 
      **must not** be mixed in a function signature. Parameters that are both 
      input and output can make the best choice unclear. Conventions consistent
      with relatd functions **must** always be followed. When adding new 
      parameters to an existing method, the established ordering convention 
      **must** be followed. Do not just stick new parameters at the end of
      the argument list.

5.5.3 Each function argument that is not a built-in type (i.e., int, double, char, etc.) **should** be passed either by reference or as a pointer to avoid unnecessary copies.

5.5.4 Each function reference or pointer argument that is not changed by the function **must** be declared "const".

5.5.6 Variable argument lists (i.e., using ellipses "...") **should** be avoided. 

      Although this is a common practice in C code, and can be done in C++ code,
      this is typically considered a dangerous carryover from C. Variadic 
      functions are not type-safe and they require tight coupling between 
      caller and callee, and can result in undefined behavior.

5.5.7 Each argument in a function declaration **must** be given a name that exactly matches the function implementation. 

      For example, use::

         void computeSomething(int op_count, int mode);

      not::

         void computeSomething(int, int);


=============================
5.6 Function implementations
=============================

5.6.1 Each function body **should** be a reasonable length to be easily understood and viewed in a text editor. Long, complex routines **should** be refactored into smaller parts when this is reasonable to increase clarity, flexibility, and the potential for code reuse.

5.6.2 Each function **should** have exactly one return point to make control logic clear.

      Functions with multiple return points tend to be a source of errors when 
      modifying code. Such routines can always be refactored to have a single 
      return point by using local scope boolean variables and/or different 
      control logic.

      A function **may** have two return points if the first return statement 
      is associated with error condition check, for example. In this case, 
      the error check **should** be performed at the start of the function body
      before other statements are reached. For example, the following is a 
      reasonable use of two function return points because the error condition
      check and the return value for successful completion are clearly visible::

         int computeSomething(int in_val)
         {
            if (in_val < 0) { return -1; }

            // ...rest of function implementation...

            return 0;
         }

5.6.3 "Sanity checks" should be performed on values of function arguments (e.g., range checking, null pointer checking, etc.) upon entry to a function. 

      This is an excellent way to provide run-time debugging capabilities in 
      code. Currently, we have a set of *assertion* macros to make syntax
      consistent. When triggered, they can emit a failed boolean expression and
      descriptive message that help to understand the violation. They are 
      active or not based on the compilation mode, either debug (active) or 
      optimized (inactive). For example::

         void doSomething(int in_val, Foo* in_foo)
         {
            ATK_ASSERT_MSG( in_val >= 0, "in_val must be positive or zero" );
            ATK_ASSERT( in_foo != NULL );

            // ...function implementation...
         }  


======================
5.7 Inline functions
======================

Function inlining is a compile time operation and the full definition of an 
inline function must be seen wherever it is called. Thus, any function to be
inlined must be implemented in a header file. 

When a function is implemented in a header file, but not declared inline, a 
compiler will choose whether or not to inline the function. Typically, 
a compiler will not inline a function that is too long or too complex (e.g.,
if it contains complicated conditional logic). When a compiler inlines a 
function, it replace the function call with the body of the function. Most
modern ccompilers do a good job of deciding when inlining is a good choice.

**Important notes:**

  * When a function implementation appears in a header file, every file that
    uses that inline method will often also emit a *function version* of the 
    method in the object file (\*.o file). This is needed to properly
    support function pointers.
  * When a function is explicitly declared inline, using the "inline" keyword,
    the compiler still decides whether to inline the function. It is possible to
    specify function attributes and compiler flags that will force a compiler to
    inline a function. Excessive inlining can cause executable code bloat and 
    may make debugging dificult. Thus, care must be used when deciding which 
    functions to explicitly declare inline. 

**When in doubt, don't use the "inline" keyword and let the compiler decide whether to inline a function.**

5.7.1 Simple, short frequently called functions, such as accessors, **should** be implemented inline in header files in most cases.

      **Exception:** Most accessors that return an object by value (i.e., not by
      pointer or a reference) **should not** be inlined. For example::

         clas MyData 
         {
            // ...public interface...
         private:
            // non-trivial private data members
            vector<Foo> m_foovec;
            Bar m_bar;
         };

         class MyClass
         {
         publis:
            //...
            MyData getData() { return m_mydata; } 

         private:
            MyData m_mydata;
         }; 

5.7.2 Class constructors **should not** be inlined. 

      A class constructor implicitly calls the constructors for its base 
      classes and initializes some or all of its data members, potentially 
      calling more constructors. If a constructor is inlined, the construction 
      and initialization needed for its members and bases will appear at every 
      object declaration.

      **Exception:** The only case where it is reasonable to inline a 
      constructor is when it has only POD ("plain old data") mambers, is not a 
      subclass of a base class, and does not explcitly declare a destructor. 
      In this case, a compiler will not even generate a destructor in most 
      cases. For example::

           class MyClass
           {
           public:
              MyClass() : m_data1(0), m_data2(0) { }

              // No destructor declared

              // ...rest of class definition...
           private:
              // class has only POD members
              int m_data1; 
              int m_data2; 
           };

5.7.3 Virtual functions **must not** be inlined due to polymorphism. 

      For example, do not declare a virtual class member function as::

         virtual void foo( ) { }

      In most circumstances, a virtual method cannot be inlined even though it
      would be inlined otherwise (e.g., because it is very short). A compiler
      must do runtime dispatch on a virtual method when it doesn't know the
      complete type at compile time.

      **Exception:** It is safe to define an empty destructor inline in an
      abstract base class with no data members. For example:: 

           class MyBase
           {
           public:
              virtual ~MyBase() {}

              virtual void doSomething(int param1) = 0;

              virtual void doSomethingElse(int param2) = 0;

              // ...

              // ...no data members...
           };


=======================================
5.8 Function and operator overloading
=======================================

5.8.1 Functions with the same name **must** differ in their argument lists and/or in their "const" attribute. 

      C++ does not allow identically named functions to differ only in their 
      return type since it is always the option of the caller to ignore or use 
      the function return value.

5.8.2 Function overloading **must not** be used to define functions that do conceptually different things. 

      Someone reading declarations of overloaded functions should be able to 
      assume (and rightfully so!) that functions with the same name do 
      something very similar.

5.8.3 If an overloaded virtual method in a base class is overridden in a derived class, all overloaded methods with the same name in the base class **must** be overridden in the derived class. 

      This prevents unexpected behavior when calling such member functions. 
      Remember that when a virtual function is overridden, the overloads of 
      that function in the base class **are not visible** to the derived class.

5.8.4 Operator overloading **must not** be used to be clever to the point of obfuscation and cause others to think too hard about an operation. Specifically, an overloaded operator must preserve "natural" semantics by appealing to common conventions and **must** have meaning similar to non-overloaded operators of the same name.

      Overloading operators can be beneficial, but **should not** be overused 
      or abused. Operator overloading is essentially "syntactic sugar" and an
      overloaded operator is just a function like any other function. An 
      important benefit of overloading is that it often allows more 
      appropriate syntax that more easily communicates the meaning of an 
      operation. The resulting code can be easier to write, maintain, and 
      understand, and it may be more efficient since it may allow the compiler
      to take advantage of longer expressions than it could otherwise.

5.8.5 Both boolean operators "==" and "!=" **should** be implemented if one of them is. 

      For consistency and correctness, the "!=" operator **should** be 
      implemented using the "==" operator implementation. For example::

         bool MyClass::operator!= (const MyClass& rhs)
         {
            return !(this == rhs);
         }

5.8.6 Standard operators, such as "&&", "||", and "," (i.e., comma), **must not** be overloaded.

      The built-in versions are treated specially by the compiler. Thus, 
      programmers cannot implement their full semantics. This can cause
      confusion. For example, the order of operand evaluation cannot be 
      guaranteed when overloading operators "&&" or "||". This may cause
      problems as someone may write code that assumes that evaluation order 
      is the same as the built-in versions.


============
5.9 Types
============

5.9.1 Behavior **should not** be selected by "switching" on the type of an 
object. 

      Good object-oriented design uses virtual functions (or templates) to 
      decide behavior. Using conditional logic (e.g., in calling code) to
      decide behavior is often unsafe and error-prone, and a clear indication 
      of poor design and improper use of the C++ type system.

5.9.2 The "bool" type **should** be used in C++ code instead of "int" for boolean true/false values.

5.9.3 The "string" type **should** be used in C++ code instead of "char\*". 

      The string type supports and optimizes many character string manipulation
      operations which can be error-prone and less efficient if implemented 
      explicitly using "char\*" and standard C library functions. Note that 
      "string" and "char\*" types are easily interchangeable, which allows C++ 
      string data to be used when interacting with C routines.

5.9.4 Class type variables **should** be defined using direct initialization instead of copy initialization to avoid unwanted and spurious type conversions and constructor calls that may be generated by compilers. 

      For example, use:: 

         std::string name("Bill");

      instead of::

         std::string name = "Bill";

      or::

         std::string name = std::string("Bill");


===================
5.10 Type casting
===================

5.10.1 C-style casts **must not** be used in C++ code. 

      All type conversions **must** be done explicitly using the named C++ 
      casting operators; i.e., "static_cast", "const_cast", "dynamic_cast", 
      "reinterpret_cast".

5.10.2 The choice to use the "static_cast" or "dynamic_cast" operator on pointers **must** consider the performance context of the code.

       The "dynamic_cast" operator is a more powerful and safer way to cast 
       pointers. However, in performance critical code, dynamic cast overhead 
       may be unacceptable. Static casts are done at compile time and are 
       essentially free at runtime whereas each dynamic cast may incur hundreds        of cycles of runtime overhead. When this choice is encountered, it may
       be wise to consider other implementation alternatives.

5.10.3 The "const_cast" operator **should** be avoided. 

       Casting away "const-ness" is often a poor programming decision and can 
       introduce errors.

       **Exception:** It may be necessary in some circumstances to cast away 
       const-ness, such as when calling const-incorrect APIs.

5.10.4 The "reinterpret_cast" **must not** be used unless absolutely necessary.

       This operator was designed to perform a low-level reinterpretation of 
       the bit pattern of an operand. This is needed only in special 
       circumstances and circumvents type safety.


================
5.11 Templates
================

5.11.1 Typically, a class (or function) template **should** be used only when the behavior of the class (or function) is completely independent of the type of the object to which it is applied. 

       Note that class member templates (e.g., member functions that are 
       templates of a class that is not a template) are often useful to 
       reduce code redundancy.

5.11.2 Generic templates that have external linkage **must** be defined in the header file where they are declared since template instantiation is a compile time operation. Thus, implementations of class templates and member templates **must** be placed in the class header file, preferably after the class definition.

5.11.3 Complete specializations of member templates or function templates **must not** appear in a header file. 

       Such methods **are not templates** and may produce link errors if their 
       definitions are seen more than once.


======================================
5.12 Conditional statements and loops
======================================

5.12.1 Curly braces **must** be used in all conditionals, loops, etc. even when the content inside the braces is a "one-liner". 

       This helps prevent coding errors and misinterpretation of intent. 
       For example, use::

          if (done) { ... }

       not::

          if (done) ...

5.12.2 One-liners **may** not be used for "if" conditionals with "else/else if"  clauses when the resulting code is clear. 

       For example, either of the following styles **may** be used::

          if (done) {
             id = 3;
          } else {
             id = 0;
          }

       or::

          if (done) { id = 3; } else { id = 0; }

5.12.3 For clarity, the shortest block of an "if/else" statement **should** come first.

5.12.4 Complex "if/else if" conditionals with many "else if" clauses **should** be avoided.

      Such statements can always be refactored using local boolean variables 
      or "switch" statements. Doing so often makes code easier to read and 
      understand and may improve performance.

5.12.5 An explicit test for zero/nonzero **must** be used in a conditional unless the tested quantity is a bool or a pointer. 

      For example, a conditional based on an integer value should use::

         if (num_lines != 0) {

      not::

         if (num_lines) {

5.12.6 A switch statement **should** use curly braces for each case and use indentation, white space, and comments for readability. Also, each case **must** contain a "break" statement and a "default" case **must** be provided to catch erroneous case values. "Fall through" cases are confusing and error-prone and so **should** be made clear in the code.

      Here is an example illustrating several preferred style practices.

.. code-block:: cpp

         switch (condition) {

            case ABC : {
               ...
               break;
            }

            case DEF :  // fall-through case
            case GHI : {
               ...
            break;
            }

            default : {
            ...
            }

         }

This code example has the following desirable properties:

   * Curly braces are used for the "switch" statement and for each case.
   * Each "case" statement is indented within the "switch" statement.
   * Blank lines are used between different cases.
   * Each case containing executable statements has a "break" statement.
   * Fall-through case is documented.
   * A "default" case is provided to catch erroneous case values.

5.12.7 The "goto" statement **should not** be used. 

      Only if alternatives are considered and determined to be less desirable, 
      should a "goto" even be contemplated.


=================
5.13 White space
=================

5.13.1 Blank lines and indentation **should** be used throughout code to enhance readability. 

      Examples of helpful white space include:

         * Between operands and operators in arithmetic expressions.
         * After reserved words, such as "while", "for", "if", "switch", etc. 
           and before the parenthesis or curly brace that follows.
         * After commas separating arguments in functions.
         * After semi-colons in for-loop expressions.
         * Before and after curly braces in almost all cases.

5.13.2 White space **must not** appear between a function name and the opening parenthesis to the argument list.  In particular, if a function call is broken across source lines, the break **must not** come between the function name and the opening parenthesis.

5.13.3 Tabs **must not** be used for indentation since this can be problematic for developers with different text editor settings.


======================
5.14 Code alignment
======================

5.14.1 Each argument in a function declaration or implementation **should** appear on its own line for clarity. The first argument **may** appear on the same line as the function name. When function areguments are placed on multiple lines, they **should** be aligned vertically for clarity.

5.14.2 All statements within a function body **should** be indented within the surrounding curly braces.

5.14.3 The start of all source lines in the same scope **should** be aligned vertically, except for continuations of previous lines.

5.14.4 If a source line is broken at a comma or semi-colon, it **must** be broken after the comma or semi-colon, not before. 

      Doing otherwise, produces code that is hard to read and can lead to 
      errors.

5.14.5 If a source line is broken at an arithmetic operator (i.e., , -, etc.), it **should** be broken after the operator, not before. 

      Doing otherwise, yields code that is harder to read and can lead to 
      errors.

5.14.6 Parentheses **should** be used in non-trivial mathematical and logical expressions to clearly indicate structure and intended order of operations and to enhance readability. 

      Do not assume everyone who looks at the code knows all the rules for 
      operator precedence.
