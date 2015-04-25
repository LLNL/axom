******************
5 Code Development
******************

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

5.4.3 Pointers **should** be set to null explicitly when memory is deallocated. Since we have not yet moved to C++11, this **must** be done using "NULL".
  
5.4.4 Data managed exclusively within C++ code **must** be allocated and deallocated using the "new" and "delete" operators. 

      The operator "new" is type-safe, simpler to use, and less error-prone 
      than the "malloc" family of C functions.  C++ new/delete operators 
      **must not** be combined with C malloc/free functions.

5.4.5 Every C++ array deallocation statement **must** include "[ ]" (i.e., "delete[ ]") to avoid memory leaks. 

      The rule of thumb is: when "[ ]" appears in the allocation, then "[ ]" 
      **must** appear in the corresponding deallocation statement.  

5.4.6 Before committing code to the source repository, one **should** use memory-checking tools to verify there are no leaks and other memory misuse.

===========================
5.5 Function declarations
===========================

Any class member function that does not change a data member of the
associated class must be declared "const".
Function arguments should be ordered the same way for all routines in a
project.
Common conventions are either to put all inputs first, then outputs, or the other way around. Input and output and outputs must not be mixed in a
function signature. Parameters that are both input and output can make the best choice unclear. Always follow conventions consistent with related
functions when possible.
When adding new parameters to an existing function, the existing ordering convention must be followed. Do not just stick new parameters at the
end.

Each function argument that is not a built-in type (i.e., int, double, char, etc.) sh
ould be passed either by reference or as a pointer to avoid unnecessary copies.
Each function reference or pointer argument that is not changed by the
function must be declared "const".
Variable argument lists (i.e., using ellipses "...") should not be used. Although
available in C++, they are generally considered a dangerous carryover from C.
Variadic functions are not type-safe, they require tight coupling between caller
and callee, and often result in undefined behavior.
Each argument in a function declaration must be given a name that matches
the function implementation. For example, use
void computeSomething(int op_count,
int mode);
not
void computeSomething(int,
int);

Function implementations
Each function body should be a reasonable length to be easily understood and
viewed in a text editor. Long, complex routines should be refactored into
smaller parts when this is reasonable to increase clarity and potentially reduce
code redundancy by promoting code reuse and flexibility.
Each function should have exactly one return point. Functions with multiple
return points tend to be a source of errors when making modifications. Such
routines can always be refactored to have a single return point by using
temporary boolean variables and/or different control logic.
A function may have two return points if the first return statement is associated with error checking,; e.g., when a state is reached beyond which it
does not make sense to continue execution of the function. In this case, error checking should be performed at the start of the function body
before other statements are reached.
For example, the following is a reasonable use of two function return points:
int computeSomething(int in_val)
{
if (in_val < 0) { return -1; }
// ...
return 0;
}

"Sanity checks" should be performed on values of function arguments (e.g.,
range checking, checking for null pointers, etc.) upon entry to a function. This
is an excellent way to provide run-time debugging capabilities in code. We need
a set of macros/functions for this.

Inline functions

Function inlining is a compile time operation; thus, inline function definitions must be seen wherever they are called (e.g., defined in a header file).
When the compiler chooses to inline a function, its code body replaces the function call where it is used. Excessive inlining may produce
unwanted executable code bloat. Thus, care must be used when deciding which functions to inline. It is important to note that compilers typically
decide whether a function is inlined. Typically, compilers will not inline functions that are too long or too complex (e.g., complicated conditional
logic).

Simple, short frequently called functions, such as accessors, should be inlined
to enhance performance.
An inline functions must be defined where it is declared in a header file. Often,
a compiler will inline small, simple functions without explicit direction from the
programmer; i.e., without using the inline keyword. To make intent clear to the
compiler, this keyword should be used.
Class constructors and destructors should not be inlined. The only case where
it is reasonable to do so is when a class has no data members and it is not a
subclass of another class. Remember that a class constructor implicitly calls
the constructors for its base classes and initializes some or all of its data
members, potentially calling more constructors. If a constructor is inlined, the
construction and initialization needed for its members and bases will appear at
every object declaration.
Virtual functions must not be inlined due to polymorphism. For example, do not
declare a class member function as
virtual void foo( ) { }

Function and operator overloading
Functions with the same name must differ in their argument lists and/or in their
"const" attribute. C++ does not allow identically named functions to differ only
in their return type since it is always the option of the caller to ignore or use the
function return value.
Function overloading must not be used to define functions that do conceptually
different things. Someone reading declarations of overloaded functions should
be able to assume (and rightfully so!) that functions with the same name do
something very similar.
If an overloaded virtual method in a base class is overridden in a derived class,
all overloaded methods with the same name in the base class must be
overridden in the derived class. This prevents unexpected behavior when
calling such member functions. Remember that when a virtual function is
overridden, the overloads of that function in the base class are not visible to
the derived class.

Operator overloading must not be used to be clever to the point of obfuscation
and cause others to think too hard about an operation. An overloaded operator
must preserve "natural" semantics by appealing to common conventions and
must have meaning similar to non-overloaded operators of the same name.
Overloading operators can be beneficial, but should not be overused or abused. Operator overloading is essentially "syntactic sugar" and an
overloaded operator is just a function like any other function. An important benefit of overloading is that it often allows more appropriate syntax
that more easily communicates the meaning of an operation. The resulting code can be easier to write, maintain, and understand, and it may be
more efficient since it may allow the compiler to take advantage of longer expressions than it could otherwise.

Both boolean operators "==" and "!=" _+should+_ be implemented if one of
them is. For consistency and correctness, the "!=" operator should be
implemented using the "==" operator implementation. For example,
bool MyClass::operator!= (const MyClass& rhs)
{
return !(this == rhs);
}
Standard operators, such as "&&", "||", and "," (i.e., comma), must not be
overloaded because the built-in versions are treated specially by the compiler.
Programmers cannot implement the full semantics of these built-in operators,
which can be confusing. For example, the order of operand evaluation cannot
be guaranteed when overloading operators "&&" or "||". This may cause
problems as someone may write code that assumes that evaluation order is the
same as the built-in versions.

Types
Behavior should not be selected by "switching" on the type of an object. Virtual
functions (or templates) should be used to decide behavior, not by conditional
logic in calling code. Doing so to customize behavior is unsafe and error-prone,
and a clear indication of poor design and improper use of the C++ type system.
The C++ "bool" type should be used instead of "int" to describe boolean
true/false values.
The C++ "string" type should be used in C++ code instead of "char*". The string
type supports and optimizes many character string manipulation operations
which can be error-prone and less efficient if implemented explicitly using

"char*" and standard C library functions. Note that "string" and "char*" types
are easily interchangeable, which allows C++ string data to be used when
interacting with C routines.
Class type variables should be defined using direct initialization instead of
copy initialization to avoid unwanted and spurious conversions and
constructor calls that may be generated by some compilers. For example, use
std::string name("Bill");
instead of
std::string name = "Bill";
or
std::string name = std::string("Bill");

Type casting
C-style casts must not be used in C++ code. All type conversions must be done
explicitly using the named C++ casting operators; i.e., "static_cast",
"const_cast", "dynamic_cast", "reinterpret_cast".
The "static_cast" operator should not be used on pointers. The "dynamic_cast"
operator is a more powerful and safer way to cast pointers. Exception: In
performance critical code, dynamic cast overhead may be unacceptable (static
casts are done at compile time and are essentially free at runtime whereas each
dynamic cast may incur hundreds of cycles of runtime overhead). In such
cases, other alternatives should be considered.
The "const_cast" operator should be avoided in almost all cases. Casting away
"const" is often a poor programming decision and can introduce errors.
However, it may be necessary in some circumstances, such as when calling
const-incorrect APIs.
The "reinterpret_cast" should be avoided unless absolutely necessary. This
operator was designed to perform a low-level reinterpretation of the bit pattern
of an operand. This is needed only in special circumstances and circumvents
type safety.

Templates

Typically, a class (or function) template should be used only when the behavior

of the class (or function) is completely independent of the type of the object to
which it is applied. Class member templates (e.g., member functions that are
templates of a class that is not a template) may be useful to reduce code
redundancy.
Generic templates that have external linkage must be defined in the header file
where they are declared since template instantiation is a compile time
operation. Implementations of class templates and member templates should b
e placed in the header file after the class definition.
Inline class methods that are templates must follow the same rules as regular
(i.e., non-template) inline methods. See Section blah…
Complete specializations of member templates or function templates must not
appear in a header file. Such methods are not templates and may produce link
errors if their definitions are seen more than once.

Conditional statements and loops
Curly braces must be used in all conditionals, loops, etc. even when the
content inside the braces is a "one-liner". This helps prevent coding errors and
misinterpretation of intent. For example, use:
if (done) { ... }
not
if (done) ...

One-liners should not be used for "if" conditionals with "else/else if" clauses.
For example, use the following form:
if (done) {
id = 3;
} else {
id = 0;
}
rather than
if (done) { id = 3; } else { id = 0; }

For clarity, the shortest block of an "if/else" statement should come first.
Complex "if/else if" conditionals with many "else if" clauses should be avoided.
Such statements can always be refactored using temporary boolean variables
or "switch" statements. Doing so makes the code easier to read and
understand and may improve performance.
An explicit test for zero/nonzero must be used in a conditional unless the
tested quantity is a bool or a pointer. For example, a conditional based on an

integer value should use:
if (num_lines != 0) {
not
if (num_lines) {

A switch statement should use curly braces for each case and use indentation
and white space for readability. Also, each case must contain a "break"
statement and a "default" case must be provided to catch erroneous case
values. "Fall through" cases are confusing and error-prone. For example:
switch (condition) {
case ABC : {
...
break;
}
case DEF :
case GHI : {
...
break;
}
default : {
...
}
}
This code example has the following desirable properties:
Curly braces are used for the "switch" statement and for the individual cases.
Each "case" statement is indented within the "switch" statement.
Blank lines are used between different cases.
Each case containing executable statements has a "break" statement.
A "default" case is provided to catch erroneous case values.

The "goto" statement should not be used. Only if alternatives are considered
and determined to be less desirable, should a "goto" even be contemplated.

White space
Blank lines and indentation should be used throughout code to enhance
readability. Examples of helpful white space include:
Between operands and operators in arithmetic expressions.
After reserved words, such as "while", "for", "if", "switch", etc. and before the parenthesis or curly brace that follows.
After commas separating arguments in functions.
After semi-colons in for-loop expressions.
Before and after curly braces in almost all cases.

White space must not appear between a function name and the opening
parenthesis to the argument list. If a function call is broken across source lines,
the break must not be between the function name and the opening parenthesis.
Tabs must not be used for indentation since this can be problematic for
developers with different text editor settings.

Code alignment
Each argument in a function declaration or implementation should appear on
its own line for clarity. The first argument may appear on the same line as the
function name, but all arguments should be aligned vertically for clarity.
All statements within a function body should be indented within the
surrounding curly braces.
The start of all source lines in the same scope should be aligned vertically,
except for continuations of previous lines.
If a source line is broken at a comma or semi-colon, it must be broken
immediately after the comma or semi-colon, not before. Doing otherwise, yields
code that is hard to read and can lead to errors.
If a source line is broken at an arithmetic operator (i.e., , -, etc.), it +should be
broken immediately after the operator, not before. Doing otherwise, yields code
that is harder to read and can lead to errors.
Parentheses should be used in non-trivial mathematical and logical
expressions to clearly indicate structure (and intended order) of operations and
to enhance readability. Do not assume everyone knows all rules for operator
precedence.


