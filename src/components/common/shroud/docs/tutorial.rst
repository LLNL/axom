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
        void TUT_function1()
        {
            Function1();
            return;
        }
    }

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

    Function1();

And the Fortran version::

    use tutorial_mod

    call function1
