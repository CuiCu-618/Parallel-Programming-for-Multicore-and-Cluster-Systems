//
// Created by mac on 01/05/2021.
//
#include "mpi.h"
#include <iomanip>
#include <iostream>

int main()
{
  double tick;
  int rank;
  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  tick = MPI_Wtick();
  if (rank == 0)
    std::cout << std::scientific << "A single MPI tick is " << tick << std::endl;
  MPI_Finalize();
  return 0;
}