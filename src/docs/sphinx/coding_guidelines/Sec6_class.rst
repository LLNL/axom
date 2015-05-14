***********************************
6 Class Design and Implementation
***********************************

===================================================
6.1 C++ class definition structure and guidelines
===================================================

This section contains guidelines for structuring a C++ class definition. 
The summary here uses numbers and text to illustrate the basic structure.
Details about individual items follow.

.. code-block:: cpp

   /* (1) Class definition preceded by documentation prologue */

   /*! 
    * \brief ...summary comment text...
    *
    * ...detailed comment text... 
    */
   class MyClass
   {

      /* (2) "friend" declarations (if needed) */

   /* (3) "public" members */
   public:

      /* (3a) static member function declarations (if needed) */

      /* (3b) public member function declarations */

   /* (4) "protected" members (rarely needed) */
   protected:
 
      /* (4a) protected member function declarations (if needed) */

   /* (5) "private members */
   private:

      /* (5a) private static data members (if needed) */

      /* (5b) private member function declarations */

      /* (5c) private data member declarations */

   };

The numbers in parentheses in the following guidelines correspond to the 
numbered items in the preceding summary.

6.1.1 Each class definition **must** be preceded by a Doxygen documentation prologue. 

      See Section 4 for details.

6.1.2 Both the opening curly brace "{" and the closing curly brace "};" for a class definition **must** be on their own source lines and must be aligned vertically with the "class" reserved word.

6.1.3 "Friend" declarations should be needed rarely if at all, but if used, they**must** appear within the body of a class definition before any class member declarations (2).

      Note that placing "friend" declarations before the "public:" keyword makes      them private, which should be the case in most circumstances. 

6.1.4 Class members **must** be declared in the following order: 

      # "public" (item 3 in summary)
      # "protected" (item 4 in summary)
      # "private" (item 5 in summary)

      That is, order members by decreasing scope of audience.

6.1.5 Static class members (methods and data) **must** be used rarely, if at all. In every case, there usage **must** be considered carefully.

      When it is determined that a static member is needed, it **must** appear 
      first in the appropriate member section. Typically, static member 
      functions **should** be "public" (item 3a in summary) and static data
      members **should** be "private" (item 5a in summary).

6.1.6 Within each set of member declarations (i.e., public, protected, private), all function declarations **must** appear before data member declarations (items 3a and 3b, 4a, 5b and 5c in summary).

6.1.7 Class data members **should** be "private" almost always. If "public" or "protected" data members are even considered, this choice **must** be reviewed carefully by other team members.

      Information hiding is an essential aspect of good software engineering 
      and private data is the best means for a class to preserve its 
      invariants. Specifically, a class should maintain control of how object 
      state can be modified to minimize side effects. In addition, restricting
      direct access to class data enforces encapsulation and facilitates 
      design changes through refactoring.

      Note that "public" and "protected" data members are not included in the 
      summary above to reinforce this guideline.

6.1.9  A class constructor that takes a single *non-default* argument, or a single argument that has a *default* value, **must** be declared "explicit". 

       This will prevent compilers from performing unexpected (and, in many
       cases, unwanted!) implicit type conversions. For example::

          class MyClass
          {
          public:
             explicit MyClass(int i, double x = 0.0);
          };

       Note that, without the explicit declaration, an implicit conversion from
       an integer to an object of type MyClass is allowed. For example::

          MyClass mc = 2;

       Clearly, this is confusing. The "explicit" keyword forces the following::

          MyClass mc(2); 

       to get the same result, which is much more clear.

6.1.10 Each class member function that does not change the object state **must** be declared "const". 

       This helps compilers detect usage errors.

6.1.11 Each class member function that returns a class data member that should not be changed by the caller **must** be declared "const" and **must** return the data member as a "const" reference or pointer.

       Often, both "const" and non-"const" versions of member access functions 
       are needed so that callers may declare the variable that holds the 
       return value with the appropriate "const-ness".

6.1.12 If a class contains nested classes or other types, these definitions **should** appear before other class members (i.e., data and functions) within the appropriate section ("public" or "private") of the enclosing class definition.

       See Section 3.8 for further guidance.

6.1.13 Each class member function and data member declaration **must** be properly documented according to the guidelines Section 4.


============================================
6.2 Class member initialization and copying
============================================

6.2.1 Every class data member **must** be initialized (using default values when appropriate) in each class constructor. That is, an initializer/initialization **must** be provided for each class data member so that every object is in a well-defined state upon construction. 

      Generally, this requires a user-defined default constructor when a class 
      has POD members. Do not assume that a compiler-generated default 
      constructor will leave any member variable in a well-defined state.

      **Exception:** A class that has no member variables, including one that 
      is derived from a base class with a default constructor that provides 
      full member initialization does not require a user-defined default 
      constructor since the compiler-generated version will suffice.

6.2.2 Data member initialization **should** be used instead of assignment in constructors, especially for small classes. 

      Initialization prevents needless run-time work and is often faster.

6.2.3 For classes with complex data members, assignment within the body of the constructor **may** be preferable. 

      If the initialization process is sufficiently complex, it **may** be 
      better to perform object initialization in a method that is called 
      after object creation, such as "init()".

6.2.4 When using initialization instead of assignment to set data member values in a constructor, the data members **should** always be initialized in the order in which they appear in the class definition. 

      Compilers adhere to this order regardless of the order that members 
      appear in the class initialization list. So you may as well agree with 
      the compiler rules and avoid potential errors when initialization of 
      one member depends on the state of another.

6.2.5 A constructor **must not** call a virtual function on any data member object since an overridden method defined in a subclass cannot be called until the object is fully constructed. 

      There is no general guarantee that data members are fully-created 
      before a constructor exits.

6.2.6 All memory allocated in a class constructor **must** be de-allocated in the class destructor. 

      The intent of constructors is to acquire resources and the intent of 
      destructors is to free those resources.  This is in the same spirit as 
      guideline 5.4.1.

6.2.7 The following guidelines for class methods that may be *automatically generated by a compiler* (i.e., default constructor, destructor, copy constructor, and copy assignment operator) **must** be followed:

      * The default constructor, copy constructor, and copy assignment operator
        **should** be declared private and unimplemented when the intent is 
        that such methods should never be called. This is a good way to help 
        a compiler to catch unintended usage. For example::

           class MyClass
           {
              // ...

           private:
              // The following methods are not implemented
              MyClass();
              MyClass(const MyClass&);
              void operator=(const MyClass&);

              // ...
           };

      * The default constructor, copy constructor, destructor, and copy 
        assignment **may** be left undeclared when the compiler-generated 
        versions are appropriate. In this case, the class header file 
        **should** contain comments indicating that the compiler-generated 
        versions of these methods will be used. 

        **Exception:** If a class inherits from a base class that declares 
        these methods private, the subclass need not declare the methods 
        private. However, a comment **should** be provided in the derived 
        class stating that the parent class enforces the non-copyable 
        properties of the class.

      * If a class is default-constructable and has POD data members, including
        raw pointers, the default constructor **must** be defined explicitly
        and the data members **must** be initialized explicitly. A 
        compiler-generated version of a default constructor will not 
        initialize such members, in general.

      * Each class **must** follow the *rule of three* which states: if the 
        destructor, copy constructor, or copy-assignment operator is 
        explicitly defined, then the others must be defined or declared 
        private and left unimplemented. In other words, compiler-generated 
        and explicit versions of these methods **must** not be mixed.

      * By convention, a functor class **should** have a copy constructor and 
        copy-assignment operator. Typically, the compiler-generated versions
        are sufficient when the class has no state or non-POD data members. 
        Since such classes are usually small and simple, the compiler-generated
        versions of these methods **may** be used without documenting the use 
        of default value semantics in the functor definition.

6.2.8 An explicit implementation of a class copy-assignment operator **must** check for assignment to self, and **must** return a reference to "\*this" after copying all data members. 

      Moreover, the *copy-and-swap* idiom **should** be used. For example::

         MyClass& MyClass::operator= (const MyClass& rhs)
         {
            if (this != &rhs) {
               MyClass temp(rhs);
               swap(temp);
            }
            return *this;
         }

      Note that the copy-and-swap idiom uses the copy constructor to create a 
      temporary object and calls a swap method provided by the class that
      copies the individual class members from the temporary object. When 
      this idiom is used, the swap method must call the swap method of each 
      of its base classes. See item 6.3.14 for additional information involving
      class inheritance.

6.2.9 A class that provides an explicit implementation of the copy constructor and copy-assignment operator **must** call the copy constructor and copy-assignment operator for each of its base classes in these methods.

================
6.3 Inheritance
================

6.3.1 Class composition **should** be used instead of inheritance to extend behavior. 

      Looser coupling between objects is typically more flexible and easier 
      to maintain and refactor.

6.3.2 Class hierarchies **should** be designed so that subclasses inherit from abstract interfaces; i.e., pure virtual base classes. 

      Inheritance is often done to reuse code that exists in a base class. 
      However, there are usually better design choices to achieve reuse. 
      Good object-oriented use of inheritance is to reuse existing *calling* 
      code by exploiting base class interfaces using polymorphism. Put another 
      way, "interface inheritance" should be used instead of "implementation 
      inheritance".

6.3.3 Deep inheritance hierarchies; i.e., more than 2 or 3 levels, **should** be avoided.

6.3.4 Multiple inheritance **should** be restricted so that only one base class contains methods that are not "pure virtual"; i.e., adhering to the Java model of inheritance is most effective for avoiding abuse of inheritance.

6.3.5 One **should not** inherit from a class that was not designed to be a base class (e.g., if it does not have a virtual destructor). 

      Doing so is bad practice and can cause problems that may not be reported 
      by a compiler; e.g., hiding base class members. To add functionality, 
      one **should** employ class composition rather than by "tweaking" an 
      existing class.

6.3.6 The destructor of a class that is designed to be a base class **must** be declared "virtual". 

      However, sometimes a destructor should not be declared virtual, such as 
      when deletion through a pointer to a base class object should be 
      disallowed.

6.3.7 "Private" and "protected" inheritance **must not** be used unless you absolutely understand the ramifications of such a choice and are sure that it will not create design and implementation problems. 

      Such a choice **must** be reviewed with team members. There almost 
      always exist better alternatives to avoid these forms of inheritance.

6.3.8 Virtual functions **should** be overridden responsibly. That is, the pre- and post-conditions, default arguments, etc. of the virtual functions should be preserved. 

      Also, the behavior of an overridden virtual function **should not** 
      deviate from the intent of the base class. Remember that derived classes 
      are subsets, not supersets, of their base classes.

6.3.9 A virtual function in a base class **should only** be defined if its behavior is always valid default behavior for *any* derived class.  

6.3.10 Inherited non-virtual methods **must not** be overloaded or hidden.

6.3.11 If a virtual function in a base class is not expected to be overridden in any derived class, then the method **should not** be declared virtual.

6.3.12 If each derived class has to provide specific behavior for a base class virtual function, then it **should** be declared *pure virtual*.

6.3.13 Virtual functions **must not** be called in a class constructor or destructor. Doing so is undefined behavior according to the C++ standard. Even if it seems to work correctly, it is fragile and potentially non-portable.

6.3.14 The copy-assignment operator of a derived class **should** use the *copy-and-swap idiom* described in item 6.2.8. 

      This is preferred to explicit assignment of base class data members 
      since the code is less coupled and since some base class data members 
      may be private. 

      Here is the preferred way to implement a swap method in a derived class::

         void DerivedClass::swap(const DerivedClass& rhs)
         {
            BaseClass::swap(rhs);
            // assign DerivedClass data members...
         }

