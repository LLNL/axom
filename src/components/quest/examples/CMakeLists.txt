#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------


make_executable( EXECUTABLE_SOURCE quest_driver.cpp
                 DEPENDS_ON 
                    quest 
                    common 
                    slic
                    slam
                 IS_EXAMPLE
                    )

make_executable( EXECUTABLE_NAME quest_sphere 
                 EXECUTABLE_SOURCE sphere.cpp 
                 DEPENDS_ON 
                    quest 
                    common 
                    slic 
                    slam
                 IS_EXAMPLE
                    )

make_executable( EXECUTABLE_NAME quest_octree_driver 
                 EXECUTABLE_SOURCE octree_driver.cpp 
                 DEPENDS_ON 
                    quest 
                    common 
                    slic 
                    slam
                 IS_EXAMPLE
                    )

if (ENABLE_MPI)

  make_executable( EXECUTABLE_NAME quest_interface 
                   EXECUTABLE_SOURCE quest_interface.cpp
                   DEPENDS_ON 
                       quest
                       common
                       slic
                   IS_EXAMPLE
                      )
  
  if (ENABLE_FORTRAN)                    
    make_executable( EXECUTABLE_NAME quest_fortran
                     EXECUTABLE_SOURCE quest_fortran.f90
                     DEPENDS_ON
                       quest
                       common
                       slic
                    IS_EXAMPLE
                       )
    set_target_properties(quest_fortran PROPERTIES LINKER_LANGUAGE Fortran)
  endif()
                      
endif()