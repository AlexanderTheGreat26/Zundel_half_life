cmake_minimum_required(VERSION 3.25)
project(Zundel_half_life)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")

add_executable(Zundel_half_life main.cpp)
