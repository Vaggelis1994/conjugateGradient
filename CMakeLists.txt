cmake_minimum_required(VERSION 3.6)
project(A8B)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp pointers_computations.cpp)
add_executable(A8B ${SOURCE_FILES})