###############
# Multimat examples
###############

#set_source_files_properties(TemplateEx.cpp PROPERTIES COMPILE_FLAGS "-S -fverbose-asm " )
#get_source_file_property(_origflags UnstructMeshField.cpp COMPILE_FLAGS)
#set_source_files_properties(UnstructMeshField.cpp "${_origflags} -fwhole-program -funroll-loops ")


set(multimat_example_headers
    helper.hpp
    )

set(multimat_example_sources
    multimat_example.cpp
    )
	
set(multimat_tests_depends_on slam multimat)

#
# Build the multimat_example library
#
blt_add_executable(
    NAME        slam_multimat_example_ex
    SOURCES     ${multimat_example_sources}
                ${multimat_example_headers}
    DEPENDS_ON  ${multimat_tests_depends_on}
    FOLDER      axom/multimat/examples
   )

   

#set(multimat_tests_extra_depends)
#if(BOOST_FOUND)
#  list(APPEND multimat_tests_extra_depends boost)
#endif()

#foreach(example_source ${example_sources})
#    get_filename_component(example_name ${example_source} NAME_WE)
#    blt_add_executable(
#        NAME        "multimat_${example_name}_ex"
#        SOURCES     ${example_source}
#        OUTPUT_DIR  ${EXAMPLE_OUTPUT_DIRECTORY}
#        DEPENDS_ON  ${multimat_tests_depends_on} ${multimat_tests_extra_depends}
#        FOLDER      axom/multimat/examples
#        )
#
#    blt_add_test(
#        NAME        multimat_${example_name}
#        COMMAND     multimat_${example_name}_ex
#        )
#endforeach()


