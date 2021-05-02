//
// Created by mac on 28/04/2021.
//

#ifndef MPI_TOOLS_H
#define MPI_TOOLS_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#define BLOCK_SIZE 2
typedef unsigned int uint;

template <class T>
void print(T x[], int N)
{
  for (uint i=0; i < N; ++i)
    std::cout << std::setprecision(3) << std::setw(6) << x[i] << " ";
  std::cout << std::endl;
}

template <class T>
void print_matrix(T A[], int N)
{
  for (uint i=0; i<N; ++i){
    for (uint j=0; j<N; ++j){
      printf("%2.2f ", A[i*N+j]);
    }
    printf("\n");
  }
}

template <class T>
void print_vectoe(T b[], int N)
{
  for (uint i=0; i<N; ++i)
      printf("%2.2f ", b[i]);
    printf("\n");
}

#endif // MPI_TOOLS_H
