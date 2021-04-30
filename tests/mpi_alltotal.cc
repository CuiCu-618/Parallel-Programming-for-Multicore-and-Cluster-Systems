//
// Created by mac on 29/04/2021.
//
#include "mpi.h"
#include "tools.h"
#include <iostream>

int main(int argv, char** argc)
{
  int rank, size;
  int *local, *global, *rbuf;
  int m;
  m = atoi(argc[1]);
  MPI_Init(&argv, &argc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  local = new int[m * size];
  global = new int[m * size];
  for (int i = 0; i < m * size; i++)
  {
    local[i] = pow(10, rank) * (i + 1);
    global[i] = 0;
  }
  for (int i = 0; i < size; i++)
  {
    if (rank == i)
      print<int>(local, size * m);
  }

  MPI_Alltoall(local, 2, MPI_INT, global, 2, MPI_INT, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  print<int>(global, size * m);
  MPI_Barrier(MPI_COMM_WORLD);

  int *scount, *sdispls, *rdispls, *rcount;
  scount = new int[size];
  rdispls = new int[size];
  sdispls = new int[size];
  rcount = new int[size];
  rbuf = new int[size * size];

  for (int i = 0; i < size * size; ++i)
  {
    rbuf[i] = 0;
  }

  for (int i = 0; i < size; ++i)
  {
    scount[i] = rank + 1;
    rcount[i] = size * scount[i];
    if (i == 0)
    {
      sdispls[i] = 0;
      rdispls[i] = 0;
    }
    else
    {
      rdispls[i] = rdispls[i - 1] + i;
      sdispls[i] = sdispls[i - 1] + rank+1;
    }
  }
  if (rank == 0)
  {
    printf("\nAlltotalv...\n");
    print<int>(scount, size);
    print<int>(sdispls, size);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Alltoallv(local, scount, sdispls, MPI_INT, rbuf, rcount, rdispls, MPI_INT, MPI_COMM_WORLD);
  print<int>(rbuf, size * size);

  delete[] local;
  delete[] global;
  delete[] rcount;
  delete[] scount;
  delete[] rdispls;
  delete[] sdispls;

  MPI_Finalize();
  return 0;
}
