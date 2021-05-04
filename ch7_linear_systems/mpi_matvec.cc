//
// Created by mac on 03/05/2021.
//
#include "mpi.h"
#include "tools.h"
#include <chrono>
#include <iostream>

void matvec_seq(const double A[], const double b[], double x[], int N);
void matvec_para(const double A[], const double b[], double x[], int N, int rank, int size);

int main(int argc, char** argv)
{
  int N, index, CYCLE, rank, size;
  double *A, *b, *x, *x_0;
  if (argc != 3)
    printf("Need args: size N, cycle CYCLE!\n");

  N = atoi(argv[1]);
  CYCLE = atoi(argv[2]);

  A = new double[N * N];
  b = new double[N];
  x = new double[N];
  x_0 = new double[N];

  memset(A, 0, N * N * sizeof(A[0]));
  memset(b, 0, N * sizeof(b[0]));
  memset(x, 0, N * sizeof(x[0]));
  memset(x_0, 0, N * sizeof(x_0[0]));

  // init A x_0
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      index = i * N + j;
      A[index] = rand() / (double)RAND_MAX + 1;
    }
    b[i] = rand() / (double)RAND_MAX + 1;
  }

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  double seq_time=0, para_time=0;
  for (int cycle = 0; cycle<CYCLE; ++cycle)
  {
    seq_time -= MPI_Wtime();
    matvec_seq(A, b, x, N);
    seq_time += MPI_Wtime();

    para_time -= MPI_Wtime();
    matvec_para(A, b, x_0, N, rank, size);
    para_time += MPI_Wtime();
  }


  if (rank == 0)
  {
    std::cout << "Compute A x b [" << CYCLE << "] times with size N = " << N << std::endl;
    printf("Seq  time: %.4f\nPare time: %.4f\n",seq_time,para_time);
    isPassed(x, x_0, N);
  }
  //  printf("Matrix A : \n");
  //  print_matrix<double>(A, N);
  //  printf("RHS b : \n");
  //  print_vectoe<double>(b, N);
  //  printf("RHS x : \n");
//    print_vectoe<double>(x_0, N);
//    print_vectoe<double>(x, N);

  delete[] A;
  delete[] b;
  delete[] x;
  delete[] x_0;

  MPI_Finalize();

  return 0;
}

void matvec_seq(const double A[], const double b[], double x[], int N)
{
  memset(x, 0, N * sizeof(x[0]));

  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      x[i] += A[i * N + j] * b[j];
    }
  }
}

void matvec_para(const double A[], const double b[], double x[], int N, int rank, int size)
{
  int row;
  /*
  double* temp;
  temp = new double[N];
  memset(temp, 0, N * sizeof(temp[0]));

  for (int i = 0; i < N / size; ++i)
  {
    row = rank + i * size;
    for (int j = 0; j < N; ++j)
    {
      temp[row] += A[row * N + j] * b[j];
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //    printf("P[%d] has row %d\n", rank, row);
  }
  for (int i = 0; i < N; i++)
  {
    MPI_Reduce(&temp[i], &x[i], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  }
   delete[] temp;

  */
  for (int i=0; i<N; ++i){
    if (i % size == rank){
      for (int j = 0; j < N; ++j)
        x[i] += A[i * N + j] * b[j];
//      MPI_Bcast(&x[i],1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
    }
  }

}