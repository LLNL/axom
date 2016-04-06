#-------------------------------------------------------------------------------
# Copyright (c) 2015, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# All rights reserved.
#
# This source code cannot be distributed without permission and further
# review from Lawrence Livermore National Laboratory.
#-------------------------------------------------------------------------------


make_executable( NAME quest_driver
                 SOURCES quest_driver.cpp
                 DEPENDS_ON 
                    quest 
                    common 
                    slic
                    slam
                 IS_EXAMPLE
                    )

make_executable( NAME quest_sphere 
                 SOURCES sphere.cpp 
                 DEPENDS_ON 
                    quest 
                    common 
                    slic 
                    slam
                 IS_EXAMPLE
                    )

make_executable( NAME quest_octree_driver 
                 SOURCES octree_driver.cpp 
                 DEPENDS_ON 
                    quest 
                    common 
                    slic 
                    slam
                 IS_EXAMPLE
                    )

if (ENABLE_MPI)

  make_executable( NAME quest_interface 
                   SOURCES quest_interface.cpp
                   DEPENDS_ON 
                       quest
                       common
                       slic
                   IS_EXAMPLE
                      )
  # Note (KW 3/2016):  There is a current problem with the fortran MPI module
  #                    on non-intel compilers (see JIRA ticket ATK-649).
  #                    For now, we are only enabling this fortran example on intel.
  if (ENABLE_FORTRAN AND COMPILER_FAMILY_IS_INTEL)                    
    make_executable( NAME quest_interface_fortran
                     SOURCES quest_interface_fortran.f90
                     DEPENDS_ON
                       quest
                       common
                       slic
                    IS_EXAMPLE
                       )
    set_target_properties(quest_interface_fortran 
                          PROPERTIES LINKER_LANGUAGE Fortran)
  endif()
                      
endif()