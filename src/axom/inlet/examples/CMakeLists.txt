# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
#------------------------------------------------------------------------------
# Examples for Inlet component
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Specify list of examples
#------------------------------------------------------------------------------

if (SOL_FOUND)
   set(inlet_examples
        arrays.cpp
        verification.cpp
        user_defined_type.cpp
        documentation_generation.cpp
        nested_structs.cpp )

    foreach(example ${inlet_examples})
        get_filename_component( ex_name ${example} NAME_WE )
        blt_add_executable( NAME       inlet_${ex_name}_ex
                            SOURCES    ${example}
                            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                            DEPENDS_ON inlet
                            FOLDER     axom/inlet/examples )
    endforeach()

    axom_add_test( NAME    inlet_verification_ex
                   COMMAND inlet_verification_ex )

    axom_add_test( NAME    inlet_documentation_generation_ex
                   COMMAND inlet_documentation_generation_ex --enableDocs --file ${CMAKE_CURRENT_LIST_DIR}/example1.lua )

    axom_add_test( NAME inlet_user_defined_type_ex
                   COMMAND inlet_user_defined_type_ex --file ${CMAKE_CURRENT_LIST_DIR}/example1.lua)

    if (MFEM_FOUND)
        blt_add_executable( NAME       inlet_mfem_coefficient_ex
                            SOURCES    mfem_coefficient.cpp
                            OUTPUT_DIR ${EXAMPLE_OUTPUT_DIRECTORY}
                            DEPENDS_ON inlet mfem
                            FOLDER     axom/inlet/examples )
        axom_add_test( NAME    inlet_mfem_coefficient_ex
                       COMMAND inlet_mfem_coefficient_ex --file ${CMAKE_CURRENT_LIST_DIR}/mfem_coef.lua)
    endif()

endif()
