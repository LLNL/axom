#include<omp.h>
#include<iostream>

/**
 * Simple test program that uses openmp without guards.  Should only be compiled when OpenMP is enabled.
 */ 
int main()
{
    int thId = omp_get_thread_num();
    int thNum = omp_get_num_threads();
    int thMax = omp_get_max_threads();
    
    std::cout <<"\nMy thread id is: " << thId 
              <<"\nNum threads is: " << thNum 
              <<"\nMax threads is: " << thMax
               << std::endl;

    return 0;
}

