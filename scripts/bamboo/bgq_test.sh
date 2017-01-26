# chang28 01-23-2017, this is the script that kicks off test
# chang28 01-25-2017, 

echo "running bgq_test.sh version 0.9"
pwd
echo "cd ../../atk_build/tests"
cd ../../atk_build/tests
echo "srun -N1 -n1 -ppsmall ./slic_asserts "
srun -N1 -n1 -ppsmall ./slic_asserts
echo "srun -N1 -n1 -ppsmall ./slic_fmt "
srun -N1 -n1 -ppsmall ./slic_fmt
