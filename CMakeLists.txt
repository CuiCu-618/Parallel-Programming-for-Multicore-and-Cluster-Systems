cmake_minimum_required(VERSION 3.11)
enable_testing()
project(MPI)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
include_directories(include)
add_subdirectory(tests)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#add_executable(MPI ex5_1.cpp)