Previous Work
=============

Communicating between language has a long history similar work.

Babel
-----

https://computation.llnl.gov/casc/components Babel parses a SIDL (Scientific Interface Definition Langauge) file to generate source. It is a hub-and-spokes approach where each language it supports is mapped to a Babel runtime object.  The last release was 2012-01-06. http://en.wikipedia.org/wiki/Babel_Middleware

Chasm
-----

http://chasm-interop.sourceforge.net/ - This page is dated July 13, 2005

Chasm is a tool to improve C++ and Fortran 90 interoperability. Chasm parses Fortran 90 source code and automatically generates C++ bridging code that can be used in C++ programs to make calls to Fortran routines. It also automatically generates C structs that provide a bridge to Fortran derived types. Chasm supplies a C++ array descriptor class which provides an interface between C and F90 arrays. This allows arrays to be created in one language and then passed to and used by the other language. http://www.cs.uoregon.edu/research/pdt/users.php

http://permalink.lanl.gov/object/tr?what=info:lanl-repo/lareport/LA-UR-01-4955

wrap
----

https://github.com/scalability-llnl/wrap

a PMPI wrapper generator

Trilinos
--------

http://trilinos.org/

Trilonos wraps C++ with C, then the Fortran over the C.  Described in the book Scientific Software Design. http://www.amazon.com/Scientific-Software-Design-The-Object-Oriented/dp/0521888131

Directory packages/ForTrilinos/src/skeleton has a basic template which must be edited to create a wrapper for a class.


Exascale Programming: Adapting What We Have Can (and Must) Work

    In 2009 and 2010, the C++ based Trilinos project developed Fortran
    interface capabilities, called ForTrilinos. As an object-oriented (OO)
    collection of libraries, we assumed that the OO features of Fortran
    2003 would provide us with natural mappings of Trilinos classes into
    Fortran equivalents. Over the two-year span of the ForTrilinos effort,
    we discovered that compiler support for 2003 features was very
    immature. ForTrilinos developers quickly came to know the handful of
    compiler developers who worked on these features and, despite close
    collaboration with them to complete and stabilize the implementation
    of Fortran 2003 features (in 2010), ForTrilinos stalled and is no
    longer developed.

http://www.hpcwire.com/2016/01/14/24151/

MPICH
-----

MPICH uses a custom perl scripts which has routine names and types in the source.

http://git.mpich.org/mpich.git/blob/HEAD:/src/binding/fortran/use_mpi/buildiface

GTK
---

gtk-fortran uses a python script which grep the C source to generate the Fortran.

https://github.com/jerryd/gtk-fortran/blob/master/src/cfwrapper.py

EOS8
----

EOS8 is written in C++ and plans to provide a Fortran binding.

    We pass around handles when possible, basically an integer id that maps to
    an actual C++ object underneath the covers. The handles are the first
    argument. There are a few places where pointers get passed around and
    reinterpret_cast'ed in the C++ layer.

    We looked at Babel but decided not to use it since it wasn't really being
    supported any longer, it was a bit heavyweight, and there was some
    run-time support needed.

    Right now the C++ headers are converted to C via a perl script with
    regexes with some hand-rolling when necessary. I tried to use clang and
    parse the headers, but never quite got it to work. The fortran is
    generated off the C wrappers. I'd like to take better advantage of native
    fortran datatypes, but aren't quite there yet.

CDI
---

CDI is a C and Fortran Interface to access Climate and NWP model Data. https://code.zmaw.de/projects/cdi

"One part of CDI[1] is a such generator. It still has some rough edges and we haven't yet decided what to do about functions returning char * (it seems like that will need some wrapping unless we simply return TYPE(c_ptr) and let the caller deal with that) but if you'd like to have a starting point in Ruby try interfaces/f2003/bindGen.rb from the tarball you can download" https://groups.google.com/d/msg/comp.lang.fortran/oadwd3HHtGA/J8DD8kGeVw8J

Links
-----

  * `Generating C Interfaces <http://fortranwiki.org/fortran/show/Generating+C+Interfaces>`_
