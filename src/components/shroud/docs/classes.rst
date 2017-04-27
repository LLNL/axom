:orphan:

.. from tutorial.rst  not sure if it deserves its own page

Classes
=======

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
constructor.  In this example the empty paren expression is required
to apply the annotation to the function instead of the result.
Likewise, ``delete`` is marked as a destructor.  These annotations
will create wrappers over the ``new`` and ``delete`` keywords.

The file ``wrapClass1.h`` will have an opaque struct for the class.
This is to allows some measure of type safety over using ``void``
pointers for every instance::

    struct s_TUT_class1;
    typedef struct s_TUT_class1 TUT_class1;


    TUT_class1 * TUT_class1_new()
    {
        Class1 *SH_rv = new Class1();
        return static_cast<TUT_class1 *>(static_cast<void *>(SH_rv));
    }

    void TUT_class1_method1(TUT_class1 * self)
    {
        Class1 *SH_this = static_cast<Class1 *>(static_cast<void *>(self));
        SH_this->Method1();
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
        rv%voidptr = c_class1_new()
    end function class1_new
    
    subroutine class1_method1(obj)
        implicit none
        class(class1) :: obj
        call c_class1_method1(obj%voidptr)
    end subroutine class1_method1


The additional C++ code to call the function::

    tutorial::Class1 *cptr = new tutorial::Class1();

    cptr->Method1();

And the Fortran version::

    type(class1) cptr

    cptr = class1_new()
    call cptr%method1
