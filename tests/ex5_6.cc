#include "mpi.h"
#include "tools.h"
#include <iostream>

struct Int_Int
{
  double val;
  int index;
};

void print_debug(Int_Int x[], size_t N)
{
  for (uint i = 0; i < N; ++i)
    std::cout << "(" << x[i].val << ", " << x[i].index << ")";
  std::cout << std::endl;
}

int main(int argc, char** argv)
{
  int rank, size;
  double local[4], global[4];
  Int_Int in[4], out[4];

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  for (int i = 0; i < 4; ++i)
  {
    local[i] = rand() / (double)RAND_MAX + pow(-rank,i+1);
    in[i].val = rand() / (double)RAND_MAX + pow(-rank,i+1);
    in[i].index = rank;
  }
//  print_debug(in, 4);
  print<double>(local,4);
//  MPI_Reduce(in, out, 4, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
  MPI_Reduce(local, global, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (rank == 0)
  {
    print<double>(global,4);
//    print_debug(out, 4);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(local,4,MPI_DOUBLE,0,MPI_COMM_WORLD);
  print<double>(local, 4);
  MPI_Finalize();
  return 0;
}