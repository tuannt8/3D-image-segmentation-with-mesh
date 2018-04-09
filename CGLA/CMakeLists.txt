cmake_minimum_required(VERSION 3.7)
project(CGLA)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES library.cpp library.h)
add_library(CGLA ${SOURCE_FILES})