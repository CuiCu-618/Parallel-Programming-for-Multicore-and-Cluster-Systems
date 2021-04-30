//
// Created by mac on 29/04/2021.
//
#include "mpi.h"
#include "tools.h"
#include <iostream>

void my_sum(double* in, double* out, int* len, MPI_Datatype* datatype)
{
  for (int i = 0; i < *len; ++i)
  {
    out[i] = in[i] + out[i];
  }
}

void maxnorm(double* in, double* out, int* len, MPI_Datatype* datatype)
{
  double in_norm = 0, out_norm = 0;
  for (int i = 0; i < *len; ++i)
  {
    in_norm += in[i] * in[i];
    out_norm += out[i] * out[i];
  }
  if (in_norm > out_norm)
  {
    for (int i = 0; i < *len; ++i)
    {
      out[i] = in[i];
    }
  }
}
int main(int argc, char** argv)
{
  int rank, size;
  double *local, *global, *my_global, *max_norm;
  double norm;
  MPI_Init(&argc, &argv);
  MPI_Op MY_SUM, MAX_NORM;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  local = new double[size];
  global = new double[size];
  my_global = new double[size];
  max_norm = new double[size];

  for (int i = 0; i < size; ++i)
  {
    local[i] = rand() / (double)RAND_MAX + pow(-rank, i + 1);
  }
  print<double>(local, size);

  MPI_Reduce(local, global, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  // create mpi op
  MPI_Op_create((MPI_User_function*)my_sum, 1, &MY_SUM);
  MPI_Op_create((MPI_User_function*)maxnorm, 1, &MAX_NORM);
  MPI_Reduce(local, my_global, size, MPI_DOUBLE, MY_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(local, max_norm, size, MPI_DOUBLE, MAX_NORM, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    print<double>(global, size);
    print<double>(my_global, size);
    print<double>(max_norm, size);
  }

  delete[] local;
  delete[] global;
  delete[] my_global;
  delete[] max_norm;
  MPI_Finalize();
  return 0;
}