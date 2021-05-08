//
// Created by mac on 08/05/2021.
//
#include "mpi.h"
#include "tools.h"
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <cmath>

int main(int argc, char** argv)
{
  int rank, size, N, it;
  double time;
  double omega, tol;
  double *A, *b, *x, *x_0;
  if (argc != 5)
  {
    printf("Input: size -N, tol -tol, omega -omega, repeat times -it\n");
    return 1;
  }

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  N = atoi(argv[1]);
  tol = atoi(argv[2]);
  omega = atof(argv[3]);
  it = atoi(argv[4]);
  if (rank == 0)
    printf("GS. Solving Ax = b with size N = %d, tol = 1e-%d\n", N, (int)tol);
  tol = pow(0.1, tol);

  A = new double[N * N];
  b = new double[N];
  x = new double[N];
  x_0 = new double[N];
  memset(x, 0, N * sizeof(x[0]));

  init_system(A, b, x_0, N);
  //  print_matrix<double>(A,N);
  double local_error, global_error;
  double local_sum, global_sum;
  int start, end, n_local = N / size, cycle;

  for (int l=0; l<it; ++l)
  {
    memset(x, 0, N * sizeof(x[0]));
    time -= MPI_Wtime();
    for (cycle = 0; cycle < MAX_CYCLE; ++cycle)
    {
      local_error = 0;
      global_error = 0;
      for (int i = 0; i < N; ++i)
      {
        local_sum = 0.0;
        start = rank * N / size;
        end = start + N / size - 1;
        for (int j = start; j <= end; ++j)
        {
          if (j != i)
            local_sum += A[i * N + j] * x[j];
        }
        MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, i / n_local, MPI_COMM_WORLD);
        if (rank == i / n_local)
        {
          // GS
          //        x[i] = 1 / (double)A[i * N + i] * (b[i] - global_sum);
          // SOR
          x[i] = omega / (double)A[i * N + i] * (b[i] - global_sum) + (1 - omega) * x[i];
          local_error = fmax(local_error, abs(x[i] - x_0[i]));
        }
      }
      MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (global_error < tol)
      {
        break;
      }
    }
    time += MPI_Wtime();
  }

  if (rank == 0)
  {
    std::cout << "Time: "  << time
              <<"\niter: " << std::setw(4) << cycle << ", error: " << std::setprecision(4)
              << std::scientific << global_error << "\n";

  }
//  print_vector<double>(x, N);
//  print_vector<double>(x_0, N);

  delete[] A;
  delete[] b;
  delete[] x;
  delete[] x_0;
  MPI_Finalize();
  return 0;
}