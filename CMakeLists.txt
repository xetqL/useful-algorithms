cmake_minimum_required(VERSION 3.17)
project(useful_algorithms)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

include_directories(include)

ENABLE_TESTING()

add_executable(serial_tests
        include/serial/algorithm.hpp
        tests/main.cpp)

