#------------------------------------------------------------------------------
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
#------------------------------------------------------------------------------

################################
# Check necessary dependencies
################################
axom_component_requires(NAME       Sidre
                        COMPONENTS SLIC
                        TPLS       Conduit)

# Note: Sidre's unordered_map is either sparsehash OR C++11's unordered_map
#       Default to std::map and warn if neither is available
if(NOT SPARSEHASH_FOUND)
  if(NOT CMAKE_CXX_STANDARD EQUAL 11)
    if(NOT CMAKE_CXX_STANDARD EQUAL 14)
      message(WARNING "No hashmap implementation available. 
        Defaulting to using an ordered map for Sidre groups.
        For improved performance, set SPARSEHASH_DIR to location of a sparsehash install
        or enable C++11 (e.g. set BLT_CXX_STD to 'c++11').")
    endif()
  endif()
endif()

################################
# Add the sidre sources
################################
set(sidre_headers
    core/Buffer.hpp
    core/Group.hpp
    core/DataStore.hpp
    core/View.hpp
    core/Attribute.hpp
    core/AttrValues.hpp
    core/Iterator.hpp
    core/MapCollection.hpp
    core/SidreTypes.hpp
    core/SidreDataTypeIds.h
    core/sidre.hpp
)

set(sidre_sources
    core/Buffer.cpp
    core/Group.cpp
    core/DataStore.cpp
    core/View.cpp
    core/Attribute.cpp
    core/AttrValues.cpp
    core/Iterator.cpp
    core/SidreUtilities.cpp
)

if(ENABLE_MPI)
   list(APPEND sidre_headers
        spio/IOManager.hpp
        spio/IOBaton.hpp
       )
   list(APPEND sidre_sources
        spio/IOManager.cpp
        spio/IOBaton.cpp
       )
endif()

if(SHROUD_FOUND)
    add_subdirectory(interface)
    add_subdirectory(spio/interface)
endif()

if (ENABLE_FORTRAN)
    # The CAPI is only built to provide a Fortran interface.  If no Fortran,
    # disable CAPI.
    list(APPEND sidre_headers
        interface/sidre.h
        interface/SidreTypes.h
        # generated headers
        interface/c_fortran/wrapSidre.h
        interface/c_fortran/wrapDataStore.h
        interface/c_fortran/wrapGroup.h
        interface/c_fortran/wrapBuffer.h
        interface/c_fortran/wrapView.h
    )
    list(APPEND sidre_sources
        # generated source
        interface/c_fortran/wrapSidre.cpp
        interface/c_fortran/wrapDataStore.cpp
        interface/c_fortran/wrapGroup.cpp
        interface/c_fortran/wrapBuffer.cpp
        interface/c_fortran/wrapView.cpp
        interface/c_fortran/wrapfsidre.F
    )

    if(ENABLE_MPI)
        list(APPEND sidre_headers
             spio/interface/c_fortran/wrapIOManager.h
        )
        list(APPEND sidre_sources
             spio/interface/c_fortran/wrapfspio.f
             spio/interface/c_fortran/wrapIOManager.cpp
        )
    endif()
endif()

################################
# Make/Install the library
################################
set(sidre_depends
    core
    conduit
    conduit_relay
    slic)

if(HDF5_FOUND)
    list(APPEND sidre_depends hdf5)
endif()

if(ENABLE_MPI)
    list(APPEND sidre_depends fmt mpi)

    if(SCR_FOUND)
        list(APPEND sidre_depends scr)
    endif()
endif()

if(SPARSEHASH_FOUND)
    list(APPEND sidre_depends sparsehash)
endif()

blt_add_library( NAME       sidre
                 SOURCES    ${sidre_sources}
                 HEADERS    ${sidre_headers}
                 DEPENDS_ON ${sidre_depends}
                 FOLDER     axom/sidre
                )

axom_install_component(NAME    sidre
                       HEADERS ${sidre_headers}
                       )
if(ENABLE_FORTRAN AND ENABLE_MPI)
    install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/axom_spio.mod DESTINATION lib/fortran)
endif()

################################
# Add examples
################################
if (AXOM_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

################################
# Add tests
################################
if (AXOM_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

################################
# Add docs
################################
if (AXOM_ENABLE_DOCS)
  add_subdirectory(docs)
endif()

################################
# Add code checks
################################
axom_add_code_checks(PREFIX   sidre
                     EXCLUDES sidre/examples/lulesh2)
