# chang28 01-23-2017, this is the script that kicks off test
# chang28 01-25-2017, 
# chang28 01-26-2017, 

echo "running bgq_test.sh version 0.95"
pwd
echo "cd ../../atk_build/tests"
cd ../../atk_build/tests
echo "srun -N1 -n1 -ppsmall ./blt_gtest_smoke"
srun -N1 -n1 -ppsmall ./blt_gtest_smoke
echo "srun -N1 -n1 -ppsmall ./common_TaskTimer"
srun -N1 -n1 -ppsmall ./common_TaskTimer
echo "srun -N1 -n1 -ppsmall ./common_fileUtilities"
srun -N1 -n1 -ppsmall ./common_fileUtilities
echo "srun -N1 -n1 -ppsmall ./common_types"
srun -N1 -n1 -ppsmall ./common_types
echo "srun -N1 -n1 -ppsmall ./common_config"
srun -N1 -n1 -ppsmall ./common_config
echo "srun -N1 -n1 -ppsmall ./compiler_flag_strict_aliasing"
srun -N1 -n1 -ppsmall ./compiler_flag_strict_aliasing
echo "srun -N1 -n1 -ppsmall ./compiler_flag_unused_var"
srun -N1 -n1 -ppsmall ./compiler_flag_unused_var
echo "srun -N1 -n1 -ppsmall ./compiler_flag_uninitialized"
srun -N1 -n1 -ppsmall ./compiler_flag_uninitialized
echo "srun -N1 -n1 -ppsmall ./compiler_flag_unused_param"
srun -N1 -n1 -ppsmall ./compiler_flag_unused_param
echo "srun -N1 -n1 -ppsmall ./compiler_flag_omp_pragma"
srun -N1 -n1 -ppsmall ./compiler_flag_omp_pragma
echo "srun -N1 -n1 -ppsmall ./conduit_smoke"
srun -N1 -n1 -ppsmall ./conduit_smoke
echo "srun -N1 -n1 -ppsmall ./fmt_smoke"
srun -N1 -n1 -ppsmall ./fmt_smoke
echo "srun -N1 -n1 -ppsmall ./mint_mesh_extent"
srun -N1 -n1 -ppsmall ./mint_mesh_extent
echo "srun -N1 -n1 -ppsmall ./quest_point"
srun -N1 -n1 -ppsmall ./quest_point
echo "srun -N1 -n1 -ppsmall ./quest_intersection"
srun -N1 -n1 -ppsmall ./quest_intersection
echo "srun -N1 -n1 -ppsmall ./quest_squared_distance"
srun -N1 -n1 -ppsmall ./quest_squared_distance
echo "srun -N1 -n1 -ppsmall ./quest_vector"
srun -N1 -n1 -ppsmall ./quest_vector
echo "srun -N1 -n1 -ppsmall ./quest_signed_distance"
srun -N1 -n1 -ppsmall ./quest_signed_distance
echo "srun -N1 -n1 -ppsmall ./quest_morton"
srun -N1 -n1 -ppsmall ./quest_morton
echo "srun -N1 -n1 -ppsmall ./quest_triangle"
srun -N1 -n1 -ppsmall ./quest_triangle
echo "srun -N1 -n1 -ppsmall ./quest_spatial_octree"
srun -N1 -n1 -ppsmall ./quest_spatial_octree
echo "srun -N1 -n1 -ppsmall ./quest_octree"
srun -N1 -n1 -ppsmall ./quest_octree
echo "srun -N1 -n1 -ppsmall ./quest_input_octree"
srun -N1 -n1 -ppsmall ./quest_inout_octree
echo "srun -N1 -n1 -ppsmall ./quest_bucket_tree"
srun -N1 -n1 -ppsmall ./quest_bucket_tree
echo "srun -N1 -n1 -ppsmall ./quest_numeric_array"
srun -N1 -n1 -ppsmall ./quest_numeric_array
echo "srun -N1 -n1 -ppsmall ./quest_boundingBox"
srun -N1 -n1 -ppsmall ./quest_boundingBox
echo "srun -N1 -n1 -ppsmall ./quest_orient"
srun -N1 -n1 -ppsmall ./quest_orient
echo "srun -N1 -n1 -ppsmall ./sidre_buffer"
srun -N1 -n1 -ppsmall ./sidre_buffer
echo "srun -N1 -n1 -ppsmall ./sidre_datastore"
srun -N1 -n1 -ppsmall ./sidre_datastore
echo "srun -N1 -n1 -ppsmall ./sidre_buffer_unit"
srun -N1 -n1 -ppsmall ./sidre_buffer_unit
echo "srun -N1 -n1 -ppsmall ./sidre_class"
srun -N1 -n1 -ppsmall ./sidre_class
echo "srun -N1 -n1 -ppsmall ./sidre_opaque"
srun -N1 -n1 -ppsmall ./sidre_opaque
echo "srun -N1 -n1 -ppsmall ./sidre_external"
srun -N1 -n1 -ppsmall ./sidre_external
echo "srun -N1 -n1 -ppsmall ./sidre_native_layout"
srun -N1 -n1 -ppsmall ./sidre_native_layout
echo "srun -N1 -n1 -ppsmall ./sidre_datastore_unit"
srun -N1 -n1 -ppsmall ./sidre_datastore_unit
echo "srun -N1 -n1 -ppsmall ./sidre_smoke"
srun -N1 -n1 -ppsmall ./sidre_smoke
echo "srun -N1 -n1 -ppsmall ./sidre_group"
srun -N1 -n1 -ppsmall ./sidre_group
echo "srun -N1 -n1 -ppsmall ./sidre_view"
srun -N1 -n1 -ppsmall ./sidre_view
echo "srun -N1 -n1 -ppsmall ./slam_ModularInt"
srun -N1 -n1 -ppsmall ./slam_ModularInt
echo "srun -N1 -n1 -ppsmall ./slam_NullSet"
srun -N1 -n1 -ppsmall ./slam_NullSet
echo "srun -N1 -n1 -ppsmall ./slam_DynamicVariableRelation"
srun -N1 -n1 -ppsmall ./slam_DynamicVariableRelation
echo "srun -N1 -n1 -ppsmall ./slam_AccessingRelationDataInMap"
srun -N1 -n1 -ppsmall ./slam_AccessingRelationDataInMap
echo "srun -N1 -n1 -ppsmall ./slam_IndirectionSet"
srun -N1 -n1 -ppsmall ./slam_IndirectionSet
echo "srun -N1 -n1 -ppsmall ./slam_Map"
srun -N1 -n1 -ppsmall ./slam_Map
echo "srun -N1 -n1 -ppsmall ./slam_utilities "
srun -N1 -n1 -ppsmall ./slam_utilities "../../src/components/slam/tests"
echo "srun -N1 -n1 -ppsmall ./slam_RangeSet"
srun -N1 -n1 -ppsmall ./slam_RangeSet
echo "srun -N1 -n1 -ppsmall ./slam_Set"
srun -N1 -n1 -ppsmall ./slam_Set
echo "srun -N1 -n1 -ppsmall ./slam_StaticConstantRelation"
srun -N1 -n1 -ppsmall ./slam_StaticConstantRelation
echo "srun -N1 -n1 -ppsmall ./slam_StaticVariableRelation"
srun -N1 -n1 -ppsmall ./slam_StaticVariableRelation
echo "srun -N1 -n1 -ppsmall ./slam_tinyHydro_unitTests"
srun -N1 -n1 -ppsmall ./slam_tinyHydro_unitTests
echo "srun -N1 -n1 -ppsmall ./slam_WindowedRangeSet"
srun -N1 -n1 -ppsmall ./slam_WindowedRangeSet
echo "srun -N1 -n1 -ppsmall ./slic_asserts "
srun -N1 -n1 -ppsmall ./slic_asserts
echo "srun -N1 -n1 -ppsmall ./slic_fmt "
srun -N1 -n1 -ppsmall ./slic_fmt
echo "cd ../examples"
cd ../examples
echo "srun -N1 -n1 -ppsmall ./slic_driver "
srun -N1 -n1 -ppsmall ./slic_driver
echo "srun -N1 -n1 -ppsmall ./slic_cpplogger_example "
srun -N1 -n1 -ppsmall ./slic_cpplogger_example
echo "srun -N1 -n1 -ppsmall ./slic_logging_example "
srun -N1 -n1 -ppsmall ./slic_logging_example
