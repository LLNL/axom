# chang28 01-23-2017, this is the script that kicks of test
echo "running bgq_test.sh"
echo "cd ../../atk_build/tests"
cd ../../atk_build/tests
echo "srun -N1 -n1 -ppsmall ./slic_fmt "
srun -N1 -n1 -ppsmall ./slic_fmt
