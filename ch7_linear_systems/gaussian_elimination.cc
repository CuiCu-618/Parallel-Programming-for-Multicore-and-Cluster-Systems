//
// Created by mac on 02/05/2021.
//
#include "tools.h"
#include <cstdio>
#include <chrono>
#include <cstdlib>
#include <iostream>

void isPassed(double x[], double x_0[], int N);
int max_col(const double A[], int k, int N);
int max_col_pointer(double **row_index, int k, int N);
void exchange_row(double A[], double b[], int r, int k, int N);
void exchange_row_pointer(double b[], int r, int k, double **row_index);

int main(int argc, char** argv)
{
  int N, index, CYCLE;
  double *A, *b, *x, *x_0, **row_index, *x_p;
  if (argc != 3)
    printf("Need args: size N, cycle CYCLE!\n");

  N = atoi(argv[1]);
  CYCLE = atoi(argv[2]);
  std::cout << "Solving Ax = b " << CYCLE << " times with size N = " << N << std::endl;

  A = new double[N * N];
  b = new double[N];
  x = new double[N];
  x_0 = new double[N];
  x_p = new double[N];
  row_index = new double*[N];

  memset(A, 0, N * N * sizeof(A[0]));
  memset(b, 0, N * sizeof(b[0]));
  memset(x, 0, N * sizeof(x[0]));
  memset(x_0, 0, N * sizeof(x_0[0]));
  memset(x_p, 0, N * sizeof(x_p[0]));
  memset(row_index, 0, N * sizeof(row_index[0]));

  // init A x_0
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      index = i * N + j;
      A[index] = rand() / (double)RAND_MAX + 1;
    }
    x_0[i] = rand() / (double)RAND_MAX + 1;
    row_index[i] = &A[i*N];
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
  //  printf("Matrix A : \n");
  //  print_matrix<double>(A, N);
  //  printf("RHS b : \n");
  //  print_vectoe<double>(b, N);
  auto start = std::chrono::steady_clock::now();
  for (int cycle = 0; cycle < CYCLE; ++cycle)
  {
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
    for (int k = N - 1; k >= 0; --k)
    {
      sum = 0.0;
      for (int i = k + 1; i < N; ++i)
      {
        sum = sum + A[k * N + i] * x[i];
      }
      x[k] = 1 / A[k * N + k] * (b[k] - sum);
    }
  }
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "\nthe actual exchange of rows\n";
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

  isPassed(x,x_0,N);


  // Pointer Forward Elimination
  start = std::chrono::steady_clock::now();
  for (int cycle = 0; cycle < CYCLE; ++cycle)
  {
    int r;
    double alpha;
    for (int k = 0; k < N - 1; ++k)
    {
      r = max_col_pointer(row_index, k, N);
      if (k != r)
        exchange_row_pointer(b, r, k, row_index);
      for (int i = k + 1; i < N; ++i)
      {
        alpha = *(row_index[i] + k) / *(row_index[k] + k);
        for (int j = k; j < N; ++j)
        {
          *(row_index[i] + j) = *(row_index[i] + j) - alpha * (*(row_index[k] + j));
        }
        b[i] = b[i] - alpha * b[k];
      }
    }
    // Backward substitution
    double sum;
    for (int k = N - 1; k >= 0; --k)
    {
      sum = 0.0;
      for (int i = k + 1; i < N; ++i)
      {
        sum = sum + *(row_index[k] + i) * x_p[i];
      }
      x_p[k] = 1 / *(row_index[k] + k) * (b[k] - sum);
    }
  }

  end = std::chrono::steady_clock::now();
  elapsed_seconds = end - start;
  std::cout << "\nusing index vectors\n";
  std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
  isPassed(x_p,x_0,N);


  //  printf("x vs x_0 : \n");
  //  print_vectoe<double>(x, N);
  //  print_vectoe<double>(x_0, N);

  delete[] A;
  delete[] b;
  delete[] x;
  delete[] x_0;
  delete[] x_p;
  delete[] row_index;
  return 0;
}

// test exchange row with pointer
/*
int main()
{
  double *A, *b, *x, **row_index;
  int N = 5, index;

  A = new double[N * N];
  b = new double[N];
  x = new double[N];
  row_index = new double*[N];
  memset(A, 0, N * N * sizeof(A[0]));
  memset(b, 0, N * sizeof(b[0]));
  memset(x, 0, N * sizeof(x[0]));
  memset(row_index, 0, N * sizeof(row_index[0]));
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      index = i * N + j;
      A[index] = rand() / (double)RAND_MAX + 1;
    }
    b[i] = rand() / (double)RAND_MAX + 1;
    row_index[i] = &A[i*N];
  }

  // Pointer Forward Elimination
  int r;
  double alpha;
  for (int k = 0; k < N - 1; ++k)
  {
    r = max_col_pointer(row_index,k,N);
    if (k != r)
      exchange_row_pointer(b,r,k,row_index);
    for (int i = k + 1; i < N; ++i)
    {
      alpha = *(row_index[i] + k) / *(row_index[k] + k);
      for (int j = k; j < N; ++j)
      {
        *(row_index[i] + j) = *(row_index[i] + j) - alpha * (*(row_index[k] + j));
      }
      b[i] = b[i] - alpha * b[k];
    }
  }
  // Backward substitution
  double sum;
  for (int k = N - 1; k >= 0; --k)
  {
    sum = 0.0;
    for (int i = k + 1; i < N; ++i)
    {
      sum = sum + *(row_index[k] + i) * x[i];
    }
    x[k] = 1/ *(row_index[k] + k) * (b[k] - sum);
  }



  printf("Matrix A : \n");
  print<double>(row_index[0],N);
  print<double>(row_index[1],N);
  print<double>(row_index[2],N);
  print<double>(row_index[3],N);
  print<double>(row_index[4],N);

  printf("RHS b : \n");
  print_vectoe<double>(b, N);
  print_vectoe<double>(x, N);


  delete[] A;
  delete[] b;
  return 0;
}
*/

void isPassed(double x[], double x_0[], int N)
{
  bool ispass = true;
  for (int i = 0; i < N; ++i)
  {
    if (abs(x[i] - x_0[i]) > 1e-10)
    {
      ispass = false;
      printf("Wrong at [%2d]\n", i);
      break;
    }
  }

  if (ispass)
    printf("Passed!\n");
  else
    printf("Failed!\n");
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

int max_col_pointer(double **row_index, int k, int N)
{
  double max = abs(*(row_index[k]+k));
  double temp;
  int r = k;
  for (int i=k+1; i<N; ++i){
    temp = abs(*(row_index[i]+k));
    if (temp > max){
      max = temp;
      r = i;
    }
  }
  return r;
}

void exchange_row_pointer(double b[], int r, int k, double **row_index)
{
  double *temp, tmp;
  temp = row_index[r];
  row_index[r] = row_index[k];
  row_index[k] = temp;

  tmp = b[r];
  b[r] = b[k];
  b[k] = tmp;
}