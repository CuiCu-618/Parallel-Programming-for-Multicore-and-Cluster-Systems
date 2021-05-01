//
// Created by mac on 01/05/2021.
//
#include "mpi.h"
#include "tools.h"
#include <cstdlib>
#include <iostream>

int main(int argc, char** argv)
{
  int rank, size;
  int my_rank, my_size;
  int* ranks;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Get the group of processes in MPI_COMM_WORLD
  MPI_Group world_group, my_group1, my_group2;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);

  ranks = new int[size / 2];
  for (int i = 0; i < size / 2; i++)
  {
    ranks[i] = 2 * i;
  }
  if (rank == 0)
    print<int>(ranks, size / 2);

  MPI_Group_size(world_group, &my_size);
  MPI_Group_rank(world_group, &my_rank);
  printf("rank / my_rank : %2d/%2d, size / my_size : %2d/%2d\n",
         rank, my_rank, size, my_size);
  MPI_Barrier(MPI_COMM_WORLD);

  // my group
  MPI_Group_incl(world_group,size/2,ranks,&my_group1);
  MPI_Group_size(my_group1, &my_size);
  MPI_Group_rank(my_group1, &my_rank);
  printf("rank / my_rank : %2d/%2d, size / my_size : %2d/%2d\n",
         rank, my_rank, size, my_size);
  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Group_excl(world_group,size/2,ranks,&my_group2);
  MPI_Group_size(my_group2, &my_size);
  MPI_Group_rank(my_group2, &my_rank);
  printf("rank / my_rank : %2d/%2d, size / my_size : %2d/%2d\n",
         rank, my_rank, size, my_size);

  //my comm
  MPI_Comm prime_comm;
  MPI_Comm_create(MPI_COMM_WORLD, my_group1, &prime_comm);
  if (MPI_COMM_NULL != prime_comm) {
    MPI_Comm_rank(prime_comm, &my_rank);
    MPI_Comm_size(prime_comm, &my_size);
  }

  printf("WORLD RANK/SIZE: %d/%d \t PRIME RANK/SIZE: %d/%d\n",
         rank, size, my_rank, my_size);

  MPI_Group_free(&world_group);
  MPI_Group_free(&my_group1);
  MPI_Group_free(&my_group2);
//  MPI_Comm_free(&prime_comm);

  MPI_Finalize();
  return 0;
}