//
// Created by mac on 28/04/2021.
//

#ifndef MPI_TOOLS_H
#define MPI_TOOLS_H

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>

#define BLOCK_SIZE 2
typedef unsigned int uint;

template <class T>
void print(T x[], size_t N)
{
  for (uint i=0; i < N; ++i)
    std::cout << std::setprecision(3) << std::setw(6) << x[i] << " ";
  std::cout << std::endl;
}

std::ostream & operator<<(std::ostream &os, float a[])
{
  for (float* p = a; *p >= 0 && *p <= 10; ++p)
    os << *p << " ";
  os << std::endl;
  return os;
}


#endif // MPI_TOOLS_H
