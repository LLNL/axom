Advanced Usage
==============

character
---------

Fortran, C, and C++ all have their own semantics for character variables.

  * Fortran ``character`` variables know their length and are blank filled
  * C ``char *`` variables are assumed to be ``NULL`` terminated.
  * C++ ``std::string`` know their own length and are ``NULL`` terminated.

It is not sufficient to pass an address between Fortran and C++ like
it is with other native types.  In order to get ideomatic behavior in
the Fortran wrappers it is often necessary to copy the values.  This
is to account for blank filled vs ``NULL`` terminated.  It is helps
support ``const`` vs non-``const`` strings.

A C 'bufferify' wrapper is created which accepts the address of the
Fortran character variable with a ``int`` argument for the declared
length of the variable (``len``) and/or a ``int`` argument for the
length with blanks trimmed off (``len_trim``).
The wrapper then uses these arguments to create a ``NULL`` terminated string
or a std::string instance.

Character Arguments
^^^^^^^^^^^^^^^^^^^

blah

Character Function
^^^^^^^^^^^^^^^^^^

blah

complex
-------


derived types
-------------

