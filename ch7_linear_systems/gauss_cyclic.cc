//
// Created by mac on 03/05/2021.
//

/// @brief Gaussian elimination with row-cyclic distribution

#include "mpi.h"
#include "tools.h"
#include <cstdlib>
#include <cstring>
#include <iostream>

int main(int argc, char** argv)
{
  int rank, size, N, index;
  double *A, *b, *x, *x_0;
  if (argc != 2)
    printf("Need arg: size N!\n");
  N = atoi(argv[1]);

  A = new double[N * N];
  b = new double[N];
  x = new double[N];
  x_0 = new double[N];
  memset(A, 0, N * N * sizeof(A[0]));
  memset(b, 0, N * sizeof(b[0]));
  memset(x, 0, N * sizeof(x[0]));
  memset(x_0, 0, N * sizeof(x_0[0]));

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (N % size != 0)
    printf("Size N should be: n * size\n");
  if (rank == 0)
    std::cout << "Solving Ax = b with size N = " << N << std::endl;

  // init A x b locally
  for (int i = 0; i < N / size; ++i)
  {
    index = rank + i * size;
    printf("Process [%d] own row : %d\n", rank, index);
    for (int j = 0; j < N; ++j)
    {
      A[index * N + j] = rand() / (double)RAND_MAX + rank + 1;
    }
  }

  if (rank==0){
    for (int i=0; i<N; ++i)
      x_0[i] = rand() / (double)RAND_MAX + rank + 1;
  }
  MPI_Bcast(x_0,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  // compute b


//  print_matrix<double>(A, N);
  print_vectoe<double>(x_0, N);

  delete[] A;
  delete[] b;
  delete[] x;
  delete[] x_0;
  MPI_Finalize();
  return 0;
}