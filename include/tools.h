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

#endif // MPI_TOOLS_H
