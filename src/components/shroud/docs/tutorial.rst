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

The **options** mapping allows the user to give information to guide
the wrapping.  **library** is used to name output files and name the
Fortran module.  **cpp_header** is the name of a C++ header file which
contains the declarations for functions to be wrapped.  **functions**
is a sequence of mappings which describe the functions to wrap.

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

To help control the scope of C names, all externals default to adding
a three letter prefix.  It defaults to the first three letters of the
**library** but may be changed by setting the option **C_prefix**.

The Fortran wrapper consists of two parts.  First is an interface
which allows Fortran to call the C routine::

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

In this case the wrapper is trivial since Fortran can call the C
routine directly.  However, it provides a place to customize the API
to coerce arguments and return values as well as other features which
will be demonstrated.

The C++ code to call the function::

    #include "tutorial.hpp"

    using namespace tutorial;
    Function1();

And the Fortran version::

    use tutorial_mod

    call function1

.. note :: rename module to just tutorial.


Arguments
---------

Integer and Real
^^^^^^^^^^^^^^^^

Integer and real types are handled using the ``iso_c_binding`` module
which match them directly to the corresponding types in C++. To wrap
``Function2``::

    double Function2(double arg1, int arg2)
    {
        return arg1 + arg2;
    }

Add the declaration to the YAML file::

    functions:
    - decl: double Function2(double arg1, int arg2)

The arguments are added to the interface for the C routine using the
``value`` attribute.  They use the ``intent(IN)`` attribute since they
are pass-by-value and cannot return a value::

        function tut_function2(arg1, arg2) result(rv) &
                bind(C, name="TUT_function2")
            use iso_c_binding
            implicit none
            real(C_DOUBLE), value, intent(IN) :: arg1
            integer(C_INT), value, intent(IN) :: arg2
            real(C_DOUBLE) :: rv
        end function tut_function2

The Fortran wrapper calls the C interface directly::

    function function2(arg1, arg2) result(rv)
        use iso_c_binding
        implicit none
        real(C_DOUBLE) :: arg1
        integer(C_INT) :: arg2
        real(C_DOUBLE) :: rv
        rv = tut_function2(arg1, arg2)
    end function function2

.. note :: add intent to wrapper

Pointer arguments
-----------------

Pointers may represent an output scalar or an array::

  - decl: int Sum(int len, int *values+dimension)


Logical
^^^^^^^

Logical variable require a conversion since they are not directly
compatible with C.  In addition, how ``.true.`` and ``.false.`` are
represented internally is compiler dependent.  So compilers use 0 for
``.false.`` while other use -1.

A simple C++ function which accepts and returns a boolean argument::

    bool Function3(bool arg)
    {
        return ! arg;
    }

Added to the YAML file as before::

    functions:
    - decl: bool Function3(bool arg)

The Fortran interface and wrapper::

        function tut_function3(arg) result(rv) &
                bind(C, name="TUT_function3")
            use iso_c_binding
            implicit none
            logical(C_BOOL), value, intent(IN) :: arg
            logical(C_BOOL) :: rv
        end function tut_function3

    function function3(arg) result(rv)
        use iso_c_binding
        implicit none
        logical :: arg
        logical :: rv
        rv = booltological(tut_function3(logicaltobool(arg)))
    end function function3

The wrapper routine uses the library function ``logicaltobool`` and
``booltological`` to use the compiler to convert between the different
kinds of logical types.  This is the first example of the wrapper
doing work to create a more idiomatic Fortran API.  It is possible to
call ``TUT_function3`` directly from Fortran, but the wrapper does the
type conversion necessary to make it easier to work within an existing
Fortran application.


Character
^^^^^^^^^

Character variables have significant differences between C and
Fortran.  The Fortran interoperabilty with C feature treat a
``character`` variable of default kind as an array of
``character(kind=C_CHAR,len=1)``.  The wrapper then deals with the C
convention of ``NULL`` termination with Fortran's blank filled.

C++ routine::

    const std::string& Function4a(
        const std::string& arg1,
        const std::string& arg2)
    {
        global_str = arg1 + arg2;
        return global_str;
    }

YAML changes::

    functions
    - decl: const std::string& Function4a(const std::string& arg1, const std::string& arg2) +pure

This is the C++ prototype with the addition of a **+pure**.  This
attribute marks the routine as Fortran ``pure`` meaning there are no
side effects.  This is necessary because the function will be called
twice.  Once to compute the length of the result and once to return
the result.

Attributes also may be added by assign new fields in **attrs**::

    - decl: const std::string& Function4a(const std::string& arg1, const std::string& arg2)
      result:
        attrs:
          pure: true

The C wrapper converts the ``std::string`` into a ``char *`` which
Fortran can deal with by assigning it to a ``type(C_PTR)``::

    const char * TUT_function4a(const char * arg1, const char * arg2)
    {
        const std::string & rv = Function4a(arg1, arg2);
        return rv.c_str();
    }

With the Fortran interface::

        pure function tut_function4a(arg1, arg2) result(rv) &
                bind(C, name="TUT_function4a")
            use iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(IN) :: arg1(*)
            character(kind=C_CHAR), intent(IN) :: arg2(*)
            type(C_PTR) rv
        end function tut_function4a

And the Fortran wrapper::

    function function4a(arg1, arg2) result(rv)
        use iso_c_binding
        implicit none
        character(*) :: arg1
        character(*) :: arg2
        character(kind=C_CHAR, len=strlen_ptr(tut_function4a(trim(arg1) // C_NULL_CHAR, trim(arg2) // C_NULL_CHAR))) :: rv
        rv = fstr(tut_function4a(  &
            trim(arg1) // C_NULL_CHAR,  &
            trim(arg2) // C_NULL_CHAR))
    end function function4a

The input arguments are trimmed of trailing blanks then concatenated
with a trailing ``NULL``.  The length of result variable ``rv`` is
computed by calling the function.  Once the result is allocated,
``tut_function4a`` is called which returns a ``type(C_PTR)``.  This
result is dereferenced by ``fstr`` and copied into ``rv``.


.. note :: create std::string from address and length?

It is possible to avoid calling the C++ function twice by passing in
another argument to hold the result.  It would be up to the caller to
ensure it is long enough.  This is done by setting the option
**F_string_result_as_arg** to true.  Like all options, it may also be
set in the global **options** and it will apply to all functions::

    - decl: const std::string& Function4b(const std::string& arg1, const std::string& arg2)
      options:
        F_string_result_as_arg: output

Only the generated wrapper is different::

    subroutine function4b(arg1, arg2, output)
        use iso_c_binding
        implicit none
        character(*) :: arg1
        character(*) :: arg2
        character(*), intent(OUT) :: output
        type(C_PTR) :: rv
        rv = tut_function4b(  &
            trim(arg1) // C_NULL_CHAR,  &
            trim(arg2) // C_NULL_CHAR)
        call FccCopyPtr(output, len(output), rv)
    end subroutine function4b

``FccCopyPtr`` is a library routine to copy the ``type(C_PTR)`` into
the character variable.

The different styles are use as::

  character(30) rv4, rv4b

  rv4 = function4a("bird", "dog")
  call function4b("bird", "dog", rv4b)



Optional Arguments
------------------

Functions with default arguments are handled by the Fortran
**optional** attribute.::

    functions:
    - decl: double Function5(double arg1 = 3.13, int arg2 = 5)

The C wrapper accepts all arguments and passes them to C++.
It is the Fortran wrapper which provides the default values, not C++.
But the end result is the same.

Fortra wrapper::

    function function5(arg1, arg2) result(rv)
        use iso_c_binding
        implicit none
        real(C_DOUBLE), optional :: arg1
        real(C_DOUBLE) :: tmp_arg1
        integer(C_INT), optional :: arg2
        integer(C_INT) :: tmp_arg2
        real(C_DOUBLE) :: rv
        if (present(arg1)) then
            tmp_arg1 = arg1
        else
            tmp_arg1 = 3.13
        endif
        if (present(arg2)) then
            tmp_arg2 = arg2
        else
            tmp_arg2 = 5
        endif
        rv = tut_function5(tmp_arg1, tmp_arg2)
    end function function5


Fortran usage::

  print *, "function5", function5()
  print *, "function5", function5(0.0d0)
  print *, "function5", function5(arg2=0)
  print *, "function5", function5(2.0d0, 2)


Overloaded Functions
--------------------

C++ allows function names to be overloaded.  Fortran supports this
using a ``generic`` interface.  The C and Fortran wrappers will
generated a wrapper for each C++ function but must mangle the name to
distinguish the names.

C++::

    void Function6(const std::string &name);
    void Function6(int indx);

By default the names are mangled by adding an index to the end. This
can be controlled by setting **function_suffix** in the YAML file::

  functions:
  - decl: void Function6(const std::string& name)
    function_suffix: _from_name
  - decl: void Function6(int indx)
    function_suffix: _from_index

The generated C wrappers uses the mangled name::

    void TUT_function6_from_name(const char * name)
    {
        Function6(name);
        return;
    }

    void TUT_function6_from_index(int indx)
    {
        Function6(indx);
        return;
    }

The generated Fortran creates routines with the same mangled name but
also creates a generic interface block to allow them to be called by
the overloaded name::

    interface function6
        module procedure function6_from_name
        module procedure function6_from_index
    end interface function6

They can be used as::

  call function6_from_name("name")
  call function6_from_index(1)
  call function6("name")
  call function6(1)



Templates
---------

C++ template are handled by creating a wrapper for each type that may
be used with the template.  The C and Fortran names are mangled by
adding a type suffix to the function name.

C++::

  template<typename ArgType>
  void Function7(ArgType arg)
  {
      return;
  }

YAML::

  - decl: void Function7(ArgType arg)
    cpp_template:
      ArgType:
        - int
        - double

C wrapper::

    void TUT_function7_int(int arg)
    {
        Function7<int>(arg);
        return;
    }
    
    void TUT_function7_double(double arg)
    {
        Function7<double>(arg);
        return;
    }

The Fortran wrapper will also generate an interface block::

    interface function7
        module procedure function7_int
        module procedure function7_double
    end interface function7


Likewise, the return type can be templated but in this case no
interface block will be generated since generic function cannot vary
only by return type.


C++::

  template<typename RetType>
  RetType Function8()
  {
      return 0;
  }

YAML::

  - decl: RetType Function8()
    cpp_template:
      RetType:
        - int
        - double

C wrapper::

    int TUT_function8_int()
    {
      int rv = Function8<int>();
      return rv;
    }

    double TUT_function8_double()
    {
      double rv = Function8<double>();
      return rv;
    }

Generic Functions
-----------------

C and C++ provide a type promotion feature when calling functions
which Fortran does not support::

    void Function9(double arg);

    Function9(1.0f);
    Function9(2.0);

When Function9 is wrapped in Fortran it may only be used with the correct arguments::

    call function9(1.)
                   1
  Error: Type mismatch in argument 'arg' at (1); passed REAL(4) to REAL(8)

It would be possible to create a version of the routine in C++ which
accepts floats, but that would require changes to the library being
wrapped.  Instead it is possible to create a generic interface to the
routine by defining which variables need their types changed.  This is
similar to templates in C++ but will only impact the Fortran wrapper.
Instead of specify the Type which changes, you specify the argument which changes::

  - decl: void Function9(double arg)
    fortran_generic:
       arg:
       -  float
       -  double

This will generate only one C wrapper which accepts a double::

  void TUT_function9(double arg)
  {
      Function9(arg);
      return;
  }

But it will generate two Fortran wrappers and a generic interface
block.  Each wrapper will coerce the argument to the correct type::

    interface function9
        module procedure function9_float
        module procedure function9_double
    end interface function9

    subroutine function9_float(arg)
        use iso_c_binding
        implicit none
        real(C_FLOAT) :: arg
        call tut_function9(real(arg, C_DOUBLE))
    end subroutine function9_float
    
    subroutine function9_double(arg)
        use iso_c_binding
        implicit none
        real(C_DOUBLE) :: arg
        call tut_function9(arg)
    end subroutine function9_double

It may now be used with single or double precision arguments::

  call function9(1.0)
  call function9(1.0d0)




Types
-----

Classes
-------

Each class is wrapped in a Fortran derived type which holds a
``type(C_PTR)`` pointer to an C++ instance of the class.  Class
methods are wrapped using Fortran's type-bound procedures.  This makes
Fortran usage very similar to C++.

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
      - decl: Class1 *new()  +constructor
      - decl: void delete()  +destructor
      - decl: void Method1()

The method ``new`` has the attribute **+constructor** to mark it as a
constructor.  It must be after the argument list to make the attribute
apply to the function as a whole instead of just the result.
Likewise, ``delete`` is marked as a destructor.

The file ``wrapClass1.h`` will have an opaque struct for the class.
This is to allows some measure of type safety over using ``void``
pointers for every instance::

    struct s_TUT_class1;
    typedef struct s_TUT_class1 TUT_class1;


    TUT_class1 * TUT_class1_new()
    {
        Class1 *selfobj = new Class1();
        return static_cast<TUT_class1 *>(static_cast<void *>(selfobj));
    }

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

