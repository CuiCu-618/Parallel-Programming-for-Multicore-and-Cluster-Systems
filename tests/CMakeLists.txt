cmake_minimum_required(VERSION 2.8)

project(MPI)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file (GLOB sources *.cc)
include_directories(../include)

foreach(ccfile ${sources})
    get_filename_component(file ${ccfile} NAME_WE)
    add_executable(${file} ${ccfile})
    add_test(${file} ${file})
endforeach()
