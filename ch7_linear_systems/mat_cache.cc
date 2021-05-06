//
// Created by mac on 06/05/2021.
//
#include "tools.h"
#include <cstdio>
#include <chrono>
#include <cstdlib>
#include <iostream>

int main(int argc, char** argv)
{
  int N, CYCLE;
  double *A, temp;
  if (argc != 3)
    printf("Need args: size N, cycle CYCLE!\n");

  N = atoi(argv[1]);
  CYCLE = atoi(argv[2]);
  std::cout << "Testing cache hit/miss with Mat read/write " << CYCLE << " times with size N = " << N << std::endl;

  A = new double[N * N];
  memset(A, 0, N * N * sizeof(A[0]));

  auto s_w_row = std::chrono::steady_clock::now();
  for (int cycle = 0; cycle < CYCLE; ++cycle)
  {
    for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        A[i * N + j] = i*j;
      }
    }
  }
  auto end_w_row = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds_w_row = end_w_row - s_w_row;

  auto s_r_row = std::chrono::steady_clock::now();
  for (int cycle = 0; cycle < CYCLE; ++cycle)
  {
    for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        temp = A[i * N + j];
      }
    }
  }
  auto end_r_row = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds_r_row = end_r_row - s_r_row;

  printf("Matrix write by row: %.2fs\nMatrix read  by row: %.2fs\n\n",
    elapsed_seconds_w_row,elapsed_seconds_r_row);

  //
  auto s_w_col = std::chrono::steady_clock::now();
  for (int cycle = 0; cycle < CYCLE; ++cycle)
  {
    for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        A[j * N + i] = i*j;
      }
    }
  }
  auto end_w_col = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds_w_col = end_w_col - s_w_col;

  auto s_r_col = std::chrono::steady_clock::now();
  for (int cycle = 0; cycle < CYCLE; ++cycle)
  {
    for (int i = 0; i < N; ++i)
    {
      for (int j = 0; j < N; ++j)
      {
        temp = A[j * N + i];
      }
    }
  }
  auto end_r_col = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds_r_col = end_r_col - s_r_col;

  printf("Matrix write by col: %.2fs\nMatrix read  by col: %.2fs\n\n",
         elapsed_seconds_w_col,elapsed_seconds_r_col);
  delete[] A;
  return 0;
}