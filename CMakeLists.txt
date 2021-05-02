cmake_minimum_required(VERSION 2.8)
enable_testing()
project(MPI)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
include_directories(include)
add_subdirectory(tests)
add_subdirectory(ch7_linear_systems)

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#add_executable(MPI ex5_1.cpp)
