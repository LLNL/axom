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

This section describes the tasks to be completed when adding a new software 
component to Axom. Apart from writing code, tasks include:

  * Setting up the appropriate directory structure
  * Modifying and adding CMake files and variables
  * Generating C and Fortran interfaces
  * Writing documentation
  * Writing tests
  * Adding a 'readme' file
  * Adding a component-specific 'uncrustify' configure file (optional)

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

To illustrate, we describe the contents of the 'sidre' component directory::

  $ cd axom/src/components/sidre
  $ ls -1
  CMakeLists.txt
  README.md
  docs
  examples
  src
  tests
  uncrustify.cfg

The names of the subdirectories are descriptive of their contents.

The 'docs' directory contains documentation for the component in 
subdirectories, such as 'doxygen' and 'sphinx' which are required. Each 
component uses doxygen to build source code documentation and sphinx to
build user documentation. Other documentation directories can be added
as needed; e.g., 'dot' for det-generated figures, 'design' for design
documents, etc.

The 'src' directory contains all header and source files for the component.
The main files (typically C++) can be organized into subdirectories
in whatever manner makes sense. As is common practice in C++ libraries,
we keep C++ header and source files in the same directories. For example, 
in sidre, these files are in a subdirectory called 'core'. 

However, **all interface files for other languages must be in a subdirectory 
called 'interface'**. For example, sidre has generated files for C, Fortran, 
python, etc. in that directory.

Each component must have a 'tests' directory that contains a comprehensive
set of unit tests.

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
'CMakeConfigureFile.cmake' that defines macro constants for third-party
library (TPL) dependencies and Axom components that help enforce consistency
for conditionally-compiled code. 

  * The name of the new component must be added to the 'COMPS' variable in that file.  
  * If a new component adds a new TPL dependency, it must be added to the 'TPL_DEPS' variable.

Then, a '#cmakedefine' definition must be added for each new macro name in the
'config.hpp.in' file in the 'axom/src/include' directory.


====================================
C and Fortran Interfaces
====================================


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


