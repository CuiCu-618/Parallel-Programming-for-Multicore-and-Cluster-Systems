//
// Created by mac on 29/04/2021.
//
#include "mpi.h"
#include "tools.h"
#include <iostream>
int main(int argc, char** argv)
{
  int rank, size;
  double *local, *global, *rbuf;
  int *rcount, *displs;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  local = new double[size];
  global = new double[size * size];
  rcount = new int[size];
  displs = new int[size];
  rbuf = new double[size * (size + 1) / 2];

  // scatter
  if (rank == 0)
  {
    printf("Scatter ... \n");
    for (int i = 0; i < size * size; ++i)
    {
      global[i] = rand() / (double)RAND_MAX + pow(-rank, i + 1);
    }
    print<double>(global, size * size);
  }
  MPI_Scatter(global, size, MPI_DOUBLE, local, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  // local compute
  print<double>(local, size);
  for (int i = 0; i < size; ++i)
  {
    local[i] = local[i] * local[i];
  }

  // Gather
  MPI_Gather(local, size, MPI_DOUBLE, global, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Gatherv
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
  {
    printf("Gatherv ...\n");
    print<double>(global, size * size);
    for (int i = 0; i < size; ++i)
    {
      rcount[i] = i + 1;
      if (i == 0)
      {
        displs[i] = 0;
      }
      else
      {
        displs[i] = displs[i - 1] + i;
      }
    }
    print<int>(rcount, size);
    print<int>(displs, size);
  }

  MPI_Gatherv(local, rank + 1, MPI_DOUBLE, rbuf, rcount, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    print<double>(rbuf, size * (size + 1) / 2);
  }

  // scatterv
  if (rank == 0)
  {
    printf("Scatterv ...\n");
    for (int i=0; i<size; ++i)
    {
      rcount[i] = size-i;
      if (i == 0)
      {
        displs[i] = 0;
      }
      else
      {
        displs[i] = displs[i - 1] + size - i + 1;
      }
    }
    print<int>(rcount, size);
    print<int>(displs, size);
  }
  MPI_Scatterv(rbuf,rcount,displs,MPI_DOUBLE,local,size-rank,MPI_DOUBLE,0,MPI_COMM_WORLD);
  print<double>(local, size-rank);

  delete[] local;
  delete[] global;
  delete[] rcount;
  delete[] displs;
  delete[] rbuf;
  MPI_Finalize();

  return 0;
}