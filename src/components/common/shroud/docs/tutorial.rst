Tutorial
========

This tutorial will walk through the steps required to create a Fortran
wrapper for a simple C++ library.

Functions
---------

The simplest item to wrap is a function in the file tutorial.hpp::

   void Function1(void);

This is wrapped using a YAML file as::

  options:
      library: Tutorial
      cpp_header: tutorial.hpp
      namespace: tutorial

  functions:
    - decl: void Function1()

.. XXX support (void)?

The **options** mapping allows the user to give information to guide the wrapping.
**library** is used to name output files and name the Fortran module.
**cpp_header** is the name of a C++ header file which contains the declarations
for functions to be wrapped.
**functions** is a sequence of mappings which describe the functions to wrap.

The generated C function in file ``wrapTutorial.cpp`` is::

    #include "wrapTutorial.h"
    #include "tutorial.hpp"

    extern "C" {
    namespace tutorial {
        void TUT_function1()
        {
            Function1();
            return;
        }
    }  // namespace tutorial
    }  // extern "C"

To help control the scope of C names, all externals default to adding a three letter prefix.
It defaults to the first three letters of the **library** but may be changed by setting 
the option **C_prefix**.

The Fortran wrapper consists of two parts.  First is an interface which allows Fortran to call the C routine::

    interface
        subroutine tut_function1() &
                bind(C, name="TUT_function1")
            use iso_c_binding
            implicit none
        end subroutine tut_function1
    end interface

The other part is a Fortran wrapper::

    subroutine function1()
        use iso_c_binding
        implicit none
        call tut_function1()
    end subroutine function1

In this case the wrapper is trivial since Fortran can call the C routine directly.  However,
it provides a place to customize the API to coerce arguments and return values as well as other
features which will be demonstrated.

The C++ code to call the function::

    #include "tutorial.hpp"

    tutorial::Function1();

And the Fortran version::

    use tutorial_mod

    call function1


Classes
-------

Each class is wrapped in a Fortran derived type which holds a ``type(C_PTR)`` pointer
to an C++ instance of the class.
Class methods are wrapped using Fortran's type-bound procedures.
This makes Fortran usage very similar to C++.

Now we'll add a simple class to the library::

    class Class1
    {
    public:
        void Method1() {};
    };

To wrap the class add the lines to the YAML file::

    classes:
      - name: Class1
        methods:
        - decl: Class1 *new+constructor
          constructor: True   # better syntax?
        - decl: void Method1()

The method ``new`` has the annotation **+constructor** to mark it as a constructor.

The file ``wrapClass1.h`` will have an opaque struct for the class.  This is to allows some
measure of type safety over using ``void`` pointers for every instance::

    #ifdef EXAMPLE_WRAPPER_IMPL
    typedef void TUT_class1;
    #else
    struct s_TUT_class1;
    typedef struct s_TUT_class1 TUT_class1;
    #endif

.. note :: When the header is used with the implementation then ``EXAMPLE_WRAPPER_IMPL`` will be defined
           and the typedef will be void.  This is simply to avoid some extra casts in the implementation.

This creates the file ``wrapClass1.cpp``::

    TUT_class1 * TUT_class1_new()
    {
        Class1 *selfobj = new Class1();
        return (TUT_class1 *) selfobj;
    }

    void TUT_class1_method1(TUT_class1 * self)
    {
        Class1 *selfobj = static_cast<Class1 *>(self);
        selfobj->Method1();
        return;
    }


    // error: invalid static_cast from type 'TUT_class1* {aka s_TUT_class1*}' to type 'tutorial::Class1*'
    // extra void * cast
    void TUT_class1_method1(TUT_class1 * self)
    {
        Class1 *selfobj = static_cast<Class1 *>(static_cast<void *>(self));
        selfobj->Method1();
        return;
    }


For Fortran a derived type is created::

    type class1
        type(C_PTR) voidptr
    contains
        procedure :: method1 => class1_method1
    end type class1

And the subroutines::

    function class1_new() result(rv)
        implicit none
        type(class1) :: rv
        rv%voidptr = tut_class1_new()
    end function class1_new
    
    subroutine class1_method1(obj)
        implicit none
        class(class1) :: obj
        call tut_class1_method1(obj%voidptr)
    end subroutine class1_method1


The additional C++ code to call the function::

    tutorial::Class1 *cptr = new tutorial::Class1();

    cptr->Method1();

And the Fortran version::

    type(class1) cptr

    cptr = class1_new()
    call cptr%method1


Arguments
---------


Types
-----

Overloaded Functions
--------------------

