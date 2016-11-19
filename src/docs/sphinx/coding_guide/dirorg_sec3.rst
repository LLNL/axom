.. ##
.. ## Copyright (c) 2016, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory.
.. ##
.. ## All rights reserved.
.. ##
.. ## This file cannot be distributed without permission and
.. ## further review from Lawrence Livermore National Laboratory.
.. ##

.. _dirorgsec-label:

=====================================
3 Directory Organization
=====================================

The goal of the guidelines in this section is to make it easy to locate a file
easily and quickly. Make it easy for your fellow developers to find stuff.

------------------------------------------
Limit scope of directory contents
------------------------------------------

3.1 The contents of each directory and file **must** be well-defined and
limited so that it can be named to clearly indicate its contents. 
The goal is to prevent directories and files from becoming bloated with 
too many or diverse concepts.


----------------------------------------------------------
Put files where it's easy to find them
----------------------------------------------------------

3.2 Header files and associated implementation files **should** reside in 
the same directory, which is a common practice for C++ libraries, unless
there is a good reason to do otherwise.

3.3 Each file **must** reside in the directory that corresponds to (and named
for) the code functionality supported by the contents of the file.

