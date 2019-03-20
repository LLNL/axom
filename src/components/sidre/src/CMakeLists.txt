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

#
# Specify all headers
#
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

#
# Specify all sources
#
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
    set(sidre_fortran_sources
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
        list(APPEND sidre_fortran_sources
             spio/interface/c_fortran/wrapfspio.f
             spio/interface/c_fortran/wrapIOManager.cpp
        )
    endif()
endif()

set(sidre_depends
    axom_utils
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


#
# make the library
#
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
                 FOLDER
                     axom/sidre
                )

install(TARGETS              sidre
        EXPORT               sidre-targets
        DESTINATION          lib
        INCLUDES DESTINATION include/sidre
        )
install(FILES       ${sidre_headers}
        DESTINATION include/sidre
        )
if(ENABLE_FORTRAN)
    install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/axom_sidre.mod DESTINATION lib/fortran)

    if(ENABLE_MPI)
        install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/axom_spio.mod DESTINATION lib/fortran)
    endif()
endif()
