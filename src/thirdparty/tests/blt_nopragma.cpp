#include<iostream>

int main()
{
    const int SZ = 100;
    int arr[SZ];

    #pragma omp parallel for
    for(int i=-0; i<SZ; ++i)
        arr[i] = i;

    std::cout <<"Value of array element 0 is " << arr[0] << std::endl;

    return 0;
}
