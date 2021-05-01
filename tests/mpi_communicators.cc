//
// Created by mac on 01/05/2021.
//
#include <iostream>
#include "mpi.h"
#include "tools.h"

int main(int argc, char **argv)
{
  int rank, size;
  int row_rank, row_size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // split
  MPI_Comm row_comm;
  int color = rank / 4;
  if (rank == 3)
    color = -1;

  MPI_Comm_split(MPI_COMM_WORLD,color,rank,&row_comm);

  MPI_Comm_rank(row_comm, &row_rank);
  MPI_Comm_size(row_comm, &row_size);
  printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d\n", rank, size, row_rank, row_size);
  MPI_Comm_free(&row_comm);

  MPI_Barrier(MPI_COMM_WORLD);



  MPI_Finalize();
  return 0;
}