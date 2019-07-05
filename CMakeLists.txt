cmake_minimum_required(VERSION 3.13)

project(Numerical_Analysis)
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(src)

add_executable(a.out main.cpp src/LinearAlgebra.cpp src/Algo.cpp)

# target_include_directories(a.out PUBLIC src)
