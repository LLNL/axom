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

.. note :: This stuff is pulled from older content we had lying around. 
           It needs to be updated, checked for correctness, gaps filled in, etc.

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


