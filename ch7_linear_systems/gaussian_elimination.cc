//
// Created by mac on 02/05/2021.
//
#include "tools.h"
#include <cstdlib>
#include <iostream>

int max_col(const double A[], int k, int N);
void exchange_row(double A[], double b[], int r, int k, int N);

int main(int argc, char** argv)
{
  int N, index;
  double *A, *b, *x, *x_0;

  N = atoi(argv[1]);
  std::cout << "Solving Ax = b with size N = " << N << std::endl;

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
    x_0[i] = rand() / (double)RAND_MAX + 1;
  }
  // compute b
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      index = i * N + j;
      b[i] += A[index] * x_0[j];
    }
  }

  printf("Matrix A : \n");
  print_matrix<double>(A, N);
  printf("RHS b : \n");
  print_vectoe<double>(b, N);

  // Forward Elimination
  int r;
  double alpha;
  for (int k = 0; k < N - 1; ++k)
  {
    r = max_col(A, k, N);
    if (k != r)
      exchange_row(A, b, r, k, N);
    for (int i = k + 1; i < N; ++i)
    {
      alpha = A[i * N + k] / A[k * N + k];
      for (int j = k; j < N; ++j)
      {
        A[i * N + j] = A[i * N + j] - alpha * A[k * N + j];
      }
      b[i] = b[i] - alpha * b[k];
    }
  }
  // Backward substitution
  double sum;
  for (int k = N - 1; k >= 0; --k) {
    sum = 0.0;
    for (int i=k+1;i<N;++i){
      sum = sum + A[k*N+i] * x[i];
    }
    x[k] = 1/A[k*N+k] * (b[k] - sum);
  }

  for (int i=0; i<N; ++i){
    if (abs(x[i]-x_0[i]) > 1e-10){
      printf("Wrong at [%2d]\n", i);
    }
  }

  printf("x vs x_0 : \n");
  print_vectoe<double>(x, N);
  print_vectoe<double>(x_0, N);

  delete[] A;
  delete[] b;
  delete[] x;
  delete[] x_0;
  return 0;
}

int max_col(const double A[], int k, int N)
{
  int index, r = k;
  double max = abs(A[k * N + k]);
  for (int i = k + 1; i < N; ++i)
  {
    index = i * N + k;
    if (abs(A[index]) > max)
    {
      max = abs(A[index]);
      r = i;
    }
  }
  return r;
}

void exchange_row(double A[], double b[], int r, int k, int N)
{
  double temp;
  for (int i = 0; i < N; ++i)
  {
    temp = A[r * N + i];
    A[r * N + i] = A[k * N + i];
    A[k * N + i] = temp;
  }
  temp = b[r];
  b[r] = b[k];
  b[k] = temp;
}