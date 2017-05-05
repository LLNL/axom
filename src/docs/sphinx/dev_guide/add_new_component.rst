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

.. _addcomponent-label:

******************************************************
Adding a New Axom Component
******************************************************

This section describes the tasks to be done when adding a new software 
component to Axom. The main tasks include:

  * Creating the appropriate directory structure
  * Modifying and adding CMake files and variables
  * Generating C and Fortran interfaces
  * Writing documentation
  * Writing tests
  * Adding a 'README' file
  * Adding a component-specific 'uncrustify' configure file (optional)

The discussion here does not contain coding guidelines. Please see
`Axom Coding Guide <../../coding_guide_docs/html/index.html>`_ for that information.

====================================
Directory Structure
====================================

In the 'axom/src/components' directory, you will find a subdirectory for
each Axom component. For example::

  $ cd axom/src/components
  $ ls -1
  CmakeLists.txt
  axom_utils
  lumberjack
  mint
  ...

All files for each component are contained in subdirectories within the
top-level component directory. 

To illustrate, we describe the contents of the *sidre* component directory::

  $ cd axom/src/components/sidre
  $ ls -1
  CMakeLists.txt
  README.md
  docs
  examples
  src
  tests
  uncrustify.cfg

The names of the subdirectories should make their contents clear.

The 'docs' directory contains documentation in subdirectories for each
type of documentation. The directories 'doxygen' and 'sphinx' are required. 
Each Axom component uses doxygen for source code documentation and sphinx 
for user documentation. Other documentation directories can be used
as needed. For example, *sidre* also contains documentation diretories: 'dot' 
for dot-generated figures, and 'design' for design documents.

The 'src' directory contains all header and source files for the component.
The main files (which are typically C++) can be organized in subdirectories
withing the 'src' directory in whatever manner makes sense. For example, in 
*sidre*, these files are in a subdirectory called 'core'. As is common 
practice for C++ libraries, we keep C++ header and source files in the same 
directories. 

**All interface files for other languages must be in a subdirectory 
called 'interface'**. To make it easy for applications written in C and
Fortran, for example, to use Axom directly in their native languages,
Axom components provide APIs in these languages. For information about
how we typically generate these APIs, see :ref:`shroudfiles-label`.

Each component must have a 'tests' directory that contains a comprehensive
set of unit tests. See :ref:`testing-label` for information about writing tests
and inserting them into our testing framework.

The 'examples' directory contains simple code examples illustrating 
component usage. A directory of examples is optional, but recommended
in most cases.

====================================
CMake Files and Variables
====================================

To properly configure and compile code for a new component, and generate 
consistent make targets, existing CMake files and variables need to be
modified in addition to adding CMake files for the new component. In this
section, we describe the changes and additions that are required.

CMake macro definitions
------------------------------

The top-level CMake directory 'axom/src/cmake' contains a file called
'CMakeConfigureFile.cmake' that defines macro constants for enabling
Axom components and setting third-party library (TPL) dependencies that 
help enforce consistency for conditionally-compiled code. When a new
component or dependency is added, that file must be modified:

  #. The name of the component must be added to the 'COMPS' variable.  
  #. If a new TPL dependency is introduced, it must be added to the 'TPL_DEPS' variable.

The CMake variables are used to generate macro contants in the Axom 
configuration header file. For each new CMake variable added, an associated
'#cmakedefine' definition must be added in the 'config.hpp.in' file in the 
'axom/src/include' directory.

Modify top-level CMakeLists.txt file
----------------------------------------

When adding a new Axom component, the file 'axom/src/components/CMakeLists.txt'
must be modified to hook the component into the CMake build configuration 
system. Specifically:

    #. Add option to enable component. For example,::

         blt_add_component(COMPONENT_NAME sidre DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})

    #. Add component dependency target by adding component name to the 'axom_components' variable.
    
Add component CMakeLists.txt files
----------------------------------------

There are several CMakeLists.txt files that must be added in various component
directories. We try to maintain consistent organization and usage across all
Axom components. To illustrate, we describe the basic contents of the 
CMakeLists.txt files in the *sidre* Axom component. See those files or those 
in other components for more details.

Top-level component directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CMakeLists.txt file in the top-level component directory, e.g., 
axom/src/components/sidre, contains the following items:

  #. Project definition; e.g.,::

       project(sidre)

  #. Checks for necessary dependencies with appropriate error or warning messages.

  #. Add subdirectories with guards as needed; e.g.,::

       add_subdirectory(src)  

     and::

       if (ENABLE_TESTS)
         add_subdirectory(tests)
       endif() 

  #. CMake exports of all component targets; e.g.,::

       install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake)

  #. Code formatting target if component-specific uncrustify configuration file
     is provided; e.g.,::

       add_code_check_targets(uncrustify.cfg) 

Component src directory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CMakeLists.txt file in the component 'src' directory defines variables for
component header files, source files, and dependencies. These CMake variable 
names have the form <component name>_<variable meaning>. So, for example,
*sidre* header file names are held in the variable 'sidre_headers'. 
The source file names are held in the variable 'sidre_sources'. Dependencies 
are held in the variable 'sidre_depends'. 

.. note:: It is important to account for all conditional inclusion of items
          in these CMake variable names. For example, a C interface is 
          generated to support a Fortran API, typically. So if Fortran is
          not enabled, it is usually not necessary to include the C header 
          files in 'sidre_headers'. Similarly, do not include items in
          the dependency variable if they are not found.

This CMakeLists.txt file also adds source subdirectories as needed 
(using the CMake 'add_subdirectory' command), adds the component as a Axom
library, and adds target definitions for dependencies. For
example, the command to add *sidre* as a library is::

  blt_add_library( NAME
                       sidre
                   SOURCES
                       "${sidre_sources}"
                       "${sidre_fortran_sources}"
                   HEADERS
                       "${sidre_headers}"
                   HEADERS_OUTPUT_SUBDIR
                       sidre
                   DEPENDS_ON
                       ${sidre_depends}
                   )

All components should follow this format.

Component examples
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CMakeLists.txt files in component examples directories define the 
following items:

  #. Variables for example source files and header files as needed
     Separate variables should be used for Fortran, C++, etc. For example,
     'example_sources' for C++, 'F_example_sources' for Fortran.

  #. An executable and test variable for each example executable to be 
     generated and each executable to be run as a test. These definitions
     use the 'blt_add_executable' and 'blt_add_test' macros, respectively.
     For example::

       blt_add_executable(NAME  <example executable name>
                          SOURCES <example source>
                          OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                          DEPENDS_ON <example dependencies>)

     and::

       blt_add_test(NAME <example executable name>
                    COMMAND <example executable name>)

     Fortran executables and tests should be guarded to prevent generation if 
     Fortran is not enabled.

Component unit tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CMakeLists.txt files in component examples directories define the 
following items:

  #. Variables for test source files as needed. Separate variables should 
     be used for Fortran, C++, etc. For example, 'gtest_sidre_tests' for
     C++ tests, 'gtest_sidre_C_tests' for C tests, and 'fruit_sidre_tests'
     for Fortran tests. Note that we use the *Google Test* framework for C
     and C++ tests and *Fruit* for Fortran tests.

  #. An executable and test variable for each test executable to be 
     generated. These variables use the 'blt_add_executable' and 
     'blt_add_test' macros, respectively, as described above.

     Fortran executables and tests should be guarded to prevent generation if 
     Fortran is not enabled.




.. _shroudfiles-label:

====================================
C and Fortran Interfaces
====================================

Typically, we use the Shroud tool to generate C and
Fortran APIs from our C++ interfaces. This makes it easy for applications 
written in those languages to use Axom directly in their native languages.
To use Shroud, create a *yaml* file in the 'interface' directory named 
For example, sidre has generated files for C, Fortran, 
python, etc. in subdirectories in the 'interface' directory.


====================================
Documentation
==================================== 

====================================
Tests
====================================


====================================
README File
====================================


====================================
Uncrustify Coniguration File
====================================



Note: when adding a new component, config.hpp.in file must be updated with 
#cmakedefine AXOM_USE_<new-component-name> 

======================================================
Adding a New Axom Component
======================================================

This section describes how to modify the Axom build system when 
adding a new component.

1. Create the folder for the new component, inside the components directory.

     `<https://lc.llnl.gov/bitbucket/projects/ATK/repos/axom/browse/src/components>`_


*  NOTES:  Create a python file to create a template for a new component.

2. Edit the **CMakeLists.txt** in the src/components directory. Use the **add_component** macro to add the new component.

      **CMakeLists.txt file:** ::

         ## add components examples

         add_component(COMPONENT_NAME common DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})
         add_component(COMPONENT_NAME slic DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})
         add_component(COMPONENT_NAME meshapi DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})
         add_component(COMPONENT_NAME sidre DEFAULT_STATE ${ENABLE_ALL_COMPONENTS})

3. Inside the **src/components/<component_name>** add a new **CMakeLists.txt** .
   Each component directory may also have **docs**, **examples**, **src** and **tests** directories.

    **Example: slic directory structure:**

.. image:: ./slic_directory.png

4. Optionally each component can have its own **uncrustify.cfg** file detailing formatting choices for the code.
   In this example, the new component Foo depends on Conduit.

    **Details of Foo's 'CMakeLists.txt:** ::


             ################################
             # Datastore
             ################################
             project(foo)


             ################################
             # Check necessary dependencies
             ################################
             if(NOT CONDUIT_FOUND)
                message(FATAL_ERROR "Foo requires Conduit. Set CONDUIT_DIR to location of built Conduit.")
             endif()


             ################################
             # Add the Foo sources
             ################################
             add_subdirectory(src)


             ################################
             # Add examples
             ################################
             if (ENABLE_EXAMPLES)
                add_subdirectory(examples)
             endif()


             ################################
             # Add tests
             ################################
             if (ENABLE_TESTS)
                add_subdirectory(tests)
             endif()

             add_code_check_targets(uncrustify.cfg)


             ################################
             # Add docs
             ################################
             if (ENABLE_DOCS)
                add_subdirectory(docs)
             endif()


             ################################
             # Create CMake importable
             # exports for all of our targets
             ################################
             install(EXPORT ${PROJECT_NAME}-targets DESTINATION lib/cmake) 

5. Create another **CMakeLists.txt** file in the *src* directory of the component.
    This contains a list of the headers, sources, and how to build them. blt_add_library
    handles building and installing the library.

    **Details of Foo's 'CMakeLists.txt:** ::

             set(foo_headers
                 Foo.hpp
                 )
             
             #
             # Specify all sources
             #
             set(foo_sources
                 Foo.cpp
                 )
             
             
             #
             # make the library
             #
             blt_add_library( NAME
                                  foo
                              SOURCES
                                  "${foo_sources}"
                              HEADERS
                     "${foo_headers}"
                              HEADERS_OUTPUT_SUBDIR
                                  foo
                              DEPENDS_ON
                                  common conduit
                              )


