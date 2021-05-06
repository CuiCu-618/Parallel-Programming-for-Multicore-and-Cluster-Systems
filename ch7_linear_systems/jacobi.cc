//
// Created by mac on 06/05/2021.
//

#include "tools.h"
#include <cstdlib>
#include <iomanip>
#include <iostream>
void jacobi_it(const double A[], const double b[], double x_0[], double x[], double x_old[], int N,
               double tol);
void GS_it(const double A[], const double b[], double x_0[], double x[], int N, double tol);
void SOR_it(const double A[], const double b[], double x_0[], double x[], int N, double tol,
            double omega);

int main(int argc, char** argv)
{
  int N;
  double omega, tol;
  double *A, *b, *x, *x_0, *x_old;
  if (argc != 4)
  {
    printf("Input: size -N, tol -tol, omega -omega\n");
    return 1;
  }
  else
  {
    N = atoi(argv[1]);
    tol = atoi(argv[2]);
    omega = atof(argv[3]);

    printf("Solving Ax = b with size N = %d, tol = 1e-%d\n", N, (int)tol);
    tol = pow(0.1, tol);
  }

  A = new double[N * N];
  b = new double[N];
  x = new double[N];
  x_0 = new double[N];
  x_old = new double[N];
  memset(x, 0, N * sizeof(x[0]));
  memset(x_old, 0, N * sizeof(x_old[0]));

  init_system(A, b, x_0, N);
//  print_matrix<double>(A,N);
  //    print_vectoe<double>(b,N);
  //    print_vectoe<double>(x_0,N);
  printf("Jacobi: \n");
  jacobi_it(A, b, x_0, x, x_old, N, tol);

  printf("GS    : \n");
  memset(x, 0, N * sizeof(x[0]));
  GS_it(A, b, x_0, x, N, tol);

  printf("SOR   : \n");
  memset(x, 0, N * sizeof(x[0]));
  SOR_it(A, b, x_0, x, N, tol,omega);

  delete[] A;
  delete[] b;
  delete[] x;
  delete[] x_0;
  return 0;
}

void jacobi_it(const double A[], const double b[], double x_0[], double x[], double x_old[], int N,
               double tol)
{
  double sum, error;
  double norm, norm_0;
  norm = norm_l2(x_0, N);

  for (int cycle = 0; cycle < MAX_CYCLE; ++cycle)
  {
    memcpy(x_old, x, N * sizeof(x[0]));
    for (int i = 0; i < N; ++i)
    {
      sum = 0.0;
      for (int j = 0; j < N; ++j)
      {
        if (j != i)
          //          sum += A[i*N+j] * x[j];
          sum += A[i * N + j] * x_old[j];
      }
      x[i] = 1 / (double)A[i * N + i] * (b[i] - sum);
    }
    norm_0 = norm_l2(x, N);
    error = abs(norm - norm_0);

    if (error < tol)
    {
      std::cout << "iter: " << std::setw(4) << cycle << ", error: " << std::setprecision(4)
                << std::scientific << error << "\n";
      break;
    }
  }
}

void GS_it(const double A[], const double b[], double x_0[], double x[], int N, double tol)
{
  double sum, error;
  double norm, norm_0;
  norm = norm_l2(x_0, N);

  for (int cycle = 0; cycle < MAX_CYCLE; ++cycle)
  {
    for (int i = 0; i < N; ++i)
    {
      sum = 0.0;
      for (int j = 0; j < N; ++j)
      {
        if (j != i)
          sum += A[i*N+j] * x[j];
      }
      x[i] = 1 / (double)A[i * N + i] * (b[i] - sum);
    }
    norm_0 = norm_l2(x, N);
    error = abs(norm - norm_0);

    if (error < tol)
    {
      std::cout << "iter: " << std::setw(4) << cycle << ", error: " << std::setprecision(4)
                << std::scientific << error << "\n";
      break;
    }
  }
}

void SOR_it(const double A[], const double b[], double x_0[], double x[], int N, double tol, double omega)
{
  double sum, error;
  double norm, norm_0;
  norm = norm_l2(x_0, N);

  for (int cycle = 0; cycle < MAX_CYCLE; ++cycle)
  {
    for (int i = 0; i < N; ++i)
    {
      sum = 0.0;
      for (int j = 0; j < N; ++j)
      {
        if (j != i)
          sum += A[i*N+j] * x[j];
      }
      x[i] = omega / (double)A[i * N + i] * (b[i] - sum) + (1-omega) * x[i];
    }
    norm_0 = norm_l2(x, N);
    error = abs(norm - norm_0);

    if (error < tol)
    {
      std::cout << "iter: " << std::setw(4) << cycle << ", error: " << std::setprecision(4)
                << std::scientific << error << ", omega = " << omega << "\n";
      break;
    }
  }
}