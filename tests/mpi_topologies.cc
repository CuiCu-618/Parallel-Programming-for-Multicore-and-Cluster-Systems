//
// Created by mac on 02/05/2021.
//
#include "mpi.h"
#include "tools.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
  int rank, size;
  int ndims, *dims, *periods, *coords, *neighbours_ranks;
  ndims = atoi(argv[1]);
  dims = new int[ndims];
  periods = new int[ndims];
  coords = new int[ndims];

  memset(dims, 0, ndims * sizeof(dims[0]));
  memset(periods, 0, ndims * sizeof(periods[0]));
  memset(coords, 0, ndims * sizeof(coords[0]));

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Dims_create(size, ndims, dims);
  if (rank == 0)
  {
    printf("Create grid of %2d processed in %dD with dims size : ", size, ndims);
    print<int>(dims, ndims);
  }

  //  virtual Cartesian grid structure
  MPI_Comm comm_2d;
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 1, &comm_2d);
  MPI_Cart_coords(comm_2d, rank, ndims, coords);
  printf("MPI rank %2d:", rank);
  print<int>(coords, ndims);
  MPI_Barrier(MPI_COMM_WORLD);

  // shift
  neighbours_ranks = new int[size];

  // Declare our neighbours
  enum DIRECTIONS {DOWN, UP, LEFT, RIGHT};
  char* neighbours_names[4] = {"down", "up", "left", "right"};
  // Let consider dims[0] = X, so the shift tells us our left and right neighbours
  MPI_Cart_shift(comm_2d, 1, 1, &neighbours_ranks[LEFT], &neighbours_ranks[RIGHT]);

  // Let consider dims[1] = Y, so the shift tells us our up and down neighbours
  MPI_Cart_shift(comm_2d, 0, 1, &neighbours_ranks[DOWN], &neighbours_ranks[UP]);

  // Get my rank in the new communicator
  int my_rank;
  MPI_Comm_rank(comm_2d, &my_rank);

  for(int i = 0; i < 4; i++)
  {
    if(neighbours_ranks[i] == MPI_PROC_NULL)
      printf("[MPI process %d] I have no %5s neighbour.\n", my_rank, neighbours_names[i]);
    else
      printf("[MPI process %d] I have  a %5s neighbour: process %d.\n", my_rank, neighbours_names[i], neighbours_ranks[i]);
  }



  delete[] dims;
  delete[] periods;
  delete[] coords;
  delete[] neighbours_ranks;

  MPI_Comm_free(&comm_2d);
  MPI_Finalize();
  return 0;
}