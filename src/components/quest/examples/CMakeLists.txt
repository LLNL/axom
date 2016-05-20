#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------

set(quest_example_depends
        common 
        slic
        slam
        quest 
        )


blt_add_executable( 
    NAME        quest_driver
    SOURCES     quest_driver.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    )

blt_add_executable( 
    NAME        quest_sphere 
    SOURCES     sphere.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    )

blt_add_executable( 
    NAME        quest_octree_driver 
    SOURCES     octree_driver.cpp 
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    )

if (ENABLE_MPI)

  blt_add_executable( 
    NAME        quest_interface 
    SOURCES     quest_interface.cpp
    OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
    DEPENDS_ON  ${quest_example_depends}
    )
    
  if (ENABLE_FORTRAN)                    
    blt_add_executable( 
        NAME        quest_interface_fortran
        SOURCES     quest_interface_fortran.f90
        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
        DEPENDS_ON  ${quest_example_depends}
        )
        
    set_target_properties(quest_interface_fortran 
                          PROPERTIES LINKER_LANGUAGE Fortran)
                          
    set_target_properties(quest_interface_fortran 
        PROPERTIES COMPILE_FLAGS ${ATK_PREPROCESS_FORTRAN})
                              
    if(MPI_Fortran_USE_MPIF)
        blt_add_target_definitions(TO quest_interface_fortran TARGET_DEFINITIONS "USE_MPIF_HEADER")
    endif()                              
                          
  endif()
                      
endif()